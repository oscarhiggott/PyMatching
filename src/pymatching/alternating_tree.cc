#include "pymatching/alternating_tree.h"
#include "pymatching/graph_fill_region.h"

#include <algorithm>
#include <iterator>
#include <utility>

pm::AltTreeEdge::AltTreeEdge() : alt_tree_node(nullptr), edge(nullptr, nullptr, 0) {
}

pm::AltTreeEdge::AltTreeEdge(AltTreeNode *alt_tree_node, const CompressedEdge &edge)
    : alt_tree_node(alt_tree_node), edge(edge) {
}

bool pm::AltTreeEdge::operator==(const pm::AltTreeEdge &rhs) const {
    return edge == rhs.edge;
}

bool pm::AltTreeEdge::operator!=(const pm::AltTreeEdge &rhs) const {
    return !(rhs == *this);
}

void pm::AltTreeNode::add_child(const pm::AltTreeEdge &child) {
    children.push_back(child);
    child.alt_tree_node->parent = {this, child.edge.reversed()};
}

pm::AltTreeNode *pm::AltTreeNode::make_child(
    pm::GraphFillRegion *child_inner_region,
    pm::GraphFillRegion *child_outer_region,
    const pm::CompressedEdge &child_inner_to_outer_edge,
    const pm::CompressedEdge &child_compressed_edge) {
    auto child = new AltTreeNode(child_inner_region, child_outer_region, child_inner_to_outer_edge);
    auto child_alt_tree_edge = AltTreeEdge(child, child_compressed_edge);
    add_child(child_alt_tree_edge);
    return child;
}

pm::AltTreeNode::AltTreeNode() : inner_region(nullptr), outer_region(nullptr), visited(false) {
}

pm::AltTreeNode::AltTreeNode(
    pm::GraphFillRegion *inner_region,
    pm::GraphFillRegion *outer_region,
    const pm::CompressedEdge &inner_to_outer_edge,
    const pm::AltTreeEdge &parent,
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

pm::AltTreeNode::AltTreeNode(
    pm::GraphFillRegion *inner_region, pm::GraphFillRegion *outer_region, const pm::CompressedEdge &inner_to_outer_edge)
    : inner_region(inner_region), outer_region(outer_region), inner_to_outer_edge(inner_to_outer_edge), visited(false) {
    inner_region->alt_tree_node = this;
    outer_region->alt_tree_node = this;
}

pm::AltTreeNode::AltTreeNode(pm::GraphFillRegion *outer_region)
    : inner_region(nullptr), outer_region(outer_region), inner_to_outer_edge(nullptr, nullptr, 0), visited(false) {
}

const pm::AltTreeNode *pm::AltTreeNode::find_root() const {
    const pm::AltTreeNode *current = this;
    while (current->parent.alt_tree_node) {
        current = current->parent.alt_tree_node;
    }
    return current;
}

pm::AltTreeNode::~AltTreeNode() {
    if (outer_region)
        outer_region->alt_tree_node = nullptr;
    if (inner_region)
        inner_region->alt_tree_node = nullptr;
}

void pm::AltTreeNode::become_root() {
    // Performs a tree rotation turning this node into the root of the tree/
    if (!parent.alt_tree_node)
        return;
    auto old_parent = parent.alt_tree_node;
    old_parent->become_root();
    parent.alt_tree_node->inner_region = inner_region;
    parent.alt_tree_node->inner_to_outer_edge = parent.edge;
    inner_region = nullptr;
    pm::unstable_erase(parent.alt_tree_node->children, [&](const pm::AltTreeEdge &x) {
        return x.alt_tree_node == this;
    });
    parent = pm::AltTreeEdge();
    add_child(pm::AltTreeEdge(old_parent, inner_to_outer_edge.reversed()));
    inner_to_outer_edge = CompressedEdge();
}

bool pm::AltTreeNode::operator==(const pm::AltTreeNode &rhs) const {
    return tree_equal(rhs);
}

bool pm::AltTreeNode::operator!=(const pm::AltTreeNode &rhs) const {
    return !(rhs == *this);
}

bool pm::AltTreeNode::tree_equal(const pm::AltTreeNode &other) const {
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

std::vector<pm::AltTreeNode *> pm::AltTreeNode::all_nodes_in_tree() {
    std::vector<pm::AltTreeNode *> all_nodes;
    size_t k = 0;
    all_nodes.push_back(this);
    while (k < all_nodes.size()) {
        for (auto &c: all_nodes[k]->children) {
            all_nodes.push_back(c.alt_tree_node);
        }
        k++;
    }
    return all_nodes;
}

pm::AltTreeNode *pm::AltTreeNode::most_recent_common_ancestor(pm::AltTreeNode &other) {
    pm::AltTreeNode *this_current = this;
    this_current->visited = true;
    pm::AltTreeNode *other_current = &other;
    other_current->visited = true;
    pm::AltTreeNode *this_parent;
    pm::AltTreeNode *other_parent;
    pm::AltTreeNode *common_ancestor;
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
    while (this_parent) {
        this_parent->visited = false;
        this_parent = this_parent->parent.alt_tree_node;
    }
    return common_ancestor;
}

pm::AltTreePruneResult::AltTreePruneResult(
    std::vector<AltTreeEdge> orphan_edges, std::vector<pm::RegionEdge> pruned_path_region_edges)
    : orphan_edges(std::move(orphan_edges)), pruned_path_region_edges(std::move(pruned_path_region_edges)) {
}

pm::AltTreePruneResult pm::AltTreeNode::prune_upward_path_stopping_before(pm::AltTreeNode *prune_parent, bool back) {
    std::vector<AltTreeEdge> orphan_edges;
    std::vector<RegionEdge> pruned_path_region_edges;
    auto current_node = this;
    if (current_node != prune_parent)
        pruned_path_region_edges.reserve(3);
    // Assumes prune_parent is an ancestor
    while (current_node != prune_parent) {
        pm::move_append(current_node->children, orphan_edges);
        if (back) {
            pruned_path_region_edges.emplace_back(current_node->inner_region, current_node->inner_to_outer_edge);
            pruned_path_region_edges.emplace_back(
                current_node->parent.alt_tree_node->outer_region, current_node->parent.edge.reversed());
        } else {
            pruned_path_region_edges.emplace_back(
                current_node->outer_region, current_node->inner_to_outer_edge.reversed());
            pruned_path_region_edges.emplace_back(current_node->inner_region, current_node->parent.edge);
        }
        pm::unstable_erase(
            current_node->parent.alt_tree_node->children, [&current_node](const AltTreeEdge &child_edge) {
                return child_edge.alt_tree_node == current_node;
            });
        current_node->outer_region->alt_tree_node = nullptr;
        current_node->inner_region->alt_tree_node = nullptr;
        auto current_alias = current_node;
        current_node = current_node->parent.alt_tree_node;
        delete current_alias;
    }
    return {orphan_edges, pruned_path_region_edges};
}
