// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "pymatching/sparse_blossom/matcher/alternating_tree.h"

#include <algorithm>
#include <iterator>
#include <utility>

#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"

using namespace pm;

AltTreeEdge::AltTreeEdge() : alt_tree_node(nullptr), edge{nullptr, nullptr, 0} {
}

AltTreeEdge::AltTreeEdge(AltTreeNode *alt_tree_node, const CompressedEdge &edge)
    : alt_tree_node(alt_tree_node), edge(edge) {
}

bool AltTreeEdge::operator==(const AltTreeEdge &rhs) const {
    return edge == rhs.edge;
}

bool AltTreeEdge::operator!=(const AltTreeEdge &rhs) const {
    return !(rhs == *this);
}

void AltTreeNode::add_child(const AltTreeEdge &child) {
    children.push_back(child);
    child.alt_tree_node->parent = {this, child.edge.reversed()};
}

AltTreeNode::AltTreeNode() : inner_region(nullptr), outer_region(nullptr), visited(false) {
}

AltTreeNode::AltTreeNode(
    GraphFillRegion *inner_region,
    GraphFillRegion *outer_region,
    const CompressedEdge &inner_to_outer_edge,
    const AltTreeEdge &parent,
    std::vector<AltTreeEdge> children)
    : inner_region(inner_region),
      outer_region(outer_region),
      inner_to_outer_edge(inner_to_outer_edge),
      parent(parent),
      children(std::move(children)),
      visited(false) {
    inner_region->alt_tree_node = this;
    outer_region->alt_tree_node = this;
}

AltTreeNode::AltTreeNode(
    GraphFillRegion *inner_region, GraphFillRegion *outer_region, const CompressedEdge &inner_to_outer_edge)
    : inner_region(inner_region), outer_region(outer_region), inner_to_outer_edge(inner_to_outer_edge), visited(false) {
    inner_region->alt_tree_node = this;
    outer_region->alt_tree_node = this;
}

AltTreeNode::AltTreeNode(GraphFillRegion *outer_region)
    : inner_region(nullptr), outer_region(outer_region), inner_to_outer_edge{nullptr, nullptr, 0}, visited(false) {
    outer_region->alt_tree_node = this;
}

const AltTreeNode *AltTreeNode::find_root() const {
    const AltTreeNode *current = this;
    while (current->parent.alt_tree_node) {
        current = current->parent.alt_tree_node;
    }
    return current;
}

void AltTreeNode::become_root() {
    // Performs a tree rotation turning this node into the root of the tree/
    if (!parent.alt_tree_node)
        return;
    auto old_parent = parent.alt_tree_node;
    old_parent->become_root();
    parent.alt_tree_node->inner_region = inner_region;
    parent.alt_tree_node->inner_to_outer_edge = parent.edge;
    inner_region = nullptr;
    unstable_erase(parent.alt_tree_node->children, [&](const AltTreeEdge &x) {
        return x.alt_tree_node == this;
    });
    parent = AltTreeEdge();
    add_child(AltTreeEdge(old_parent, inner_to_outer_edge.reversed()));
    inner_to_outer_edge = CompressedEdge{nullptr, nullptr, 0};
}

bool AltTreeNode::operator==(const AltTreeNode &rhs) const {
    return tree_equal(rhs);
}

bool AltTreeNode::operator!=(const AltTreeNode &rhs) const {
    return !(rhs == *this);
}

bool AltTreeNode::tree_equal(const AltTreeNode &other) const {
    if (inner_region != other.inner_region || outer_region != other.outer_region ||
        inner_to_outer_edge != other.inner_to_outer_edge || children.size() != other.children.size()) {
        return false;
    }
    if (children.empty())
        return true;
    for (size_t i = 0; i < children.size(); i++) {
        if (children[i].edge != other.children[i].edge)
            return false;
        if (!children[i].alt_tree_node->tree_equal(*other.children[i].alt_tree_node))
            return false;
    }
    return true;
}

std::vector<AltTreeNode *> AltTreeNode::all_nodes_in_tree() {
    std::vector<AltTreeNode *> all_nodes;
    size_t k = 0;
    all_nodes.push_back(this);
    while (k < all_nodes.size()) {
        for (auto &c : all_nodes[k]->children) {
            all_nodes.push_back(c.alt_tree_node);
        }
        k++;
    }
    return all_nodes;
}

AltTreeNode *AltTreeNode::most_recent_common_ancestor(AltTreeNode &other) {
    AltTreeNode *this_current = this;
    this_current->visited = true;
    AltTreeNode *other_current = &other;
    other_current->visited = true;
    AltTreeNode *this_parent;
    AltTreeNode *other_parent;
    AltTreeNode *common_ancestor;
    while (true) {
        this_parent = this_current->parent.alt_tree_node;
        other_parent = other_current->parent.alt_tree_node;
        if (this_parent || other_parent) {
            if (this_parent) {
                this_current = this_parent;
                if (this_current->visited) {
                    common_ancestor = this_current;
                    break;
                }
                this_current->visited = true;
            }
            if (other_parent) {
                other_current = other_parent;
                if (other_current->visited) {
                    common_ancestor = other_current;
                    break;
                }
                other_current->visited = true;
            }
        } else {
            return nullptr;
        }
    }
    // Clean up 'visited' flags for the common ancestor, and its ancestors
    common_ancestor->visited = false;
    this_parent = common_ancestor->parent.alt_tree_node;
    while (this_parent && this_parent->visited) {
        this_parent->visited = false;
        this_parent = this_parent->parent.alt_tree_node;
    }
    return common_ancestor;
}

AltTreePruneResult::AltTreePruneResult(
    std::vector<AltTreeEdge> orphan_edges, std::vector<RegionEdge> pruned_path_region_edges)
    : orphan_edges(std::move(orphan_edges)), pruned_path_region_edges(std::move(pruned_path_region_edges)) {
}

AltTreePruneResult AltTreeNode::prune_upward_path_stopping_before(
    Arena<AltTreeNode> &arena, AltTreeNode *prune_parent, bool back) {
    std::vector<AltTreeEdge> orphan_edges;
    std::vector<RegionEdge> pruned_path_region_edges;
    auto current_node = this;
    if (current_node != prune_parent)
        pruned_path_region_edges.reserve(3);
    // Assumes prune_parent is an ancestor
    while (current_node != prune_parent) {
        move_append(current_node->children, orphan_edges);
        if (back) {
            pruned_path_region_edges.push_back({current_node->inner_region, current_node->inner_to_outer_edge});
            pruned_path_region_edges.push_back(
                {current_node->parent.alt_tree_node->outer_region, current_node->parent.edge.reversed()});
        } else {
            pruned_path_region_edges.push_back(
                {current_node->outer_region, current_node->inner_to_outer_edge.reversed()});
            pruned_path_region_edges.push_back({current_node->inner_region, current_node->parent.edge});
        }
        unstable_erase(current_node->parent.alt_tree_node->children, [&current_node](const AltTreeEdge &child_edge) {
            return child_edge.alt_tree_node == current_node;
        });
        current_node->outer_region->alt_tree_node = nullptr;
        current_node->inner_region->alt_tree_node = nullptr;
        auto to_remove = current_node;
        current_node = current_node->parent.alt_tree_node;
        arena.del(to_remove);
    }
    return {orphan_edges, pruned_path_region_edges};
}
