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

#ifndef PYMATCHING2_ALTERNATING_TREE_H
#define PYMATCHING2_ALTERNATING_TREE_H

#include <vector>

#include "pymatching/sparse_blossom/arena.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/region_edge.h"

namespace pm {

class AltTreeNode;
class GraphFillRegion;

struct AltTreeEdge {
    AltTreeNode* alt_tree_node;
    pm::CompressedEdge edge;
    AltTreeEdge();

    bool operator==(const AltTreeEdge& rhs) const;

    bool operator!=(const AltTreeEdge& rhs) const;

    AltTreeEdge(AltTreeNode* alt_tree_node, const CompressedEdge& edge);
};

template <class T, class UnaryPredicate>
bool unstable_erase(std::vector<T>& vec, UnaryPredicate pred);

template <class T, class UnaryPredicate>
inline bool unstable_erase(std::vector<T>& vec, UnaryPredicate pred) {
    auto res = std::find_if(vec.begin(), vec.end(), pred);
    if (res == vec.end())
        return false;
    if (vec.size() > 1) {
        *res = std::move(vec.back());
        vec.pop_back();
        return true;
    }
    vec.clear();
    return true;
}

template <typename T>
void move_append(std::vector<T>& src, std::vector<T>& dst) {
    if (dst.empty()) {
        dst = std::move(src);
    } else {
        dst.reserve(dst.size() + src.size());
        std::move(std::begin(src), std::end(src), std::back_inserter(dst));
        src.clear();
    }
}

struct AltTreePruneResult {
    std::vector<AltTreeEdge> orphan_edges;
    std::vector<RegionEdge> pruned_path_region_edges;

    AltTreePruneResult(std::vector<AltTreeEdge> orphan_edges, std::vector<pm::RegionEdge> pruned_path_region_edges);
};

/// An alternating tree is a tree with 2-colored nodes where one color grows and the other shrinks.
///
/// Additionally, an alternating tree requires that the leaves and the root are growing.
/// Additionally, it is required that each shrinking node has exactly 1 child.
///
/// These additional constraints allow an optimization used here, where each AltTreeNode is actually
/// two nodes from the alternating tree: a growing node and its single shrinking parent (the parent
/// is null for the root node).
class AltTreeNode {
   public:
    /// The shrinking node in this combined alternating tree double-node.
    GraphFillRegion* inner_region;
    /// The growing node in this combined alternating tree double-node.
    GraphFillRegion* outer_region;
    /// The edge from the shrinking region to the growing region.
    CompressedEdge inner_to_outer_edge;
    /// The edge from the shrinking region to its parent (i.e. from this double node to its parent).
    AltTreeEdge parent;
    /// The children of the growing node (i.e. the children of this double node).
    std::vector<AltTreeEdge> children;
    /// Ephemeral state used during algorithms.
    bool visited;

    AltTreeNode();
    AltTreeNode(
        GraphFillRegion* inner_region,
        GraphFillRegion* outer_region,
        const CompressedEdge& inner_to_outer_edge,
        const AltTreeEdge& parent,
        std::vector<AltTreeEdge> children);
    AltTreeNode(
        pm::GraphFillRegion* inner_region,
        pm::GraphFillRegion* outer_region,
        const pm::CompressedEdge& inner_to_outer_edge);
    AltTreeNode(pm::GraphFillRegion* outer_region);

    /// Two trees are equal if they have the same structure and they refer to the same underlying
    /// regions and have the same compressed edges.
    bool operator==(const AltTreeNode& rhs) const;
    bool operator!=(const AltTreeNode& rhs) const;

    /// Perform a tree rotation where the growing node of this double node becomes the root
    /// of the tree. The undirected edges in the underlying alternating tree do not change; only
    /// the direction of the edges (all leading away from the root) is changed.
    void become_root();
    /// Finds the most recent common ancestor between this node and the other node.
    ///
    /// Would be const except it uses ephemeral state to go faster.
    /// More specifically, breadcrumbs are left (by setting `visited=true`) as the
    /// tree is traversed upward from `this` and `other`. These are
    /// reset (`visited=false`) before the method completes for the
    /// common ancestor and its ancestors, but *not* reset for nodes
    /// on the paths from `this` and `other` stopping *before* ancestor.
    /// This does not affect the algorithm, as these nodes are always removed
    /// from the tree immediately after the common ancestor is found (either
    /// a blossom is formed, or the two trees shatter into matches).
    AltTreeNode* most_recent_common_ancestor(AltTreeNode& other);
    void add_child(const AltTreeEdge& child);
    AltTreePruneResult prune_upward_path_stopping_before(
        Arena<AltTreeNode>& arena, AltTreeNode* prune_parent, bool back);
    const AltTreeNode* find_root() const;

    /// Helper method for operator==.
    bool tree_equal(const pm::AltTreeNode& other) const;

    /// Lists all the nodes in the tree. Not used during the algorithm; might be useful for tests
    /// at some point.
    std::vector<pm::AltTreeNode*> all_nodes_in_tree();
};

}  // namespace pm

#endif  // PYMATCHING2_ALTERNATING_TREE_H
