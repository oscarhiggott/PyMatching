#ifndef PYMATCHING2_ALTERNATING_TREE_H
#define PYMATCHING2_ALTERNATING_TREE_H

#include "compressed_edge.h"
#include "graph_fill_region.h"
#include "region_edge.h"

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

class AltTreeNode {
   public:
    GraphFillRegion* inner_region;
    GraphFillRegion* outer_region;
    CompressedEdge inner_to_outer_edge;
    AltTreeEdge parent;
    std::vector<AltTreeEdge> children;  // Maybe make linked list?
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
    ~AltTreeNode();

    bool operator==(const AltTreeNode& rhs) const;

    bool operator!=(const AltTreeNode& rhs) const;

    void become_root();
    std::vector<GraphFillRegion> shatter_into_matches();
    AltTreeNode* most_recent_common_ancestor(AltTreeNode& other);
    void add_child(const AltTreeEdge& child);
    AltTreeNode* make_child(
        GraphFillRegion* child_inner_region,
        GraphFillRegion* child_outer_region,
        const CompressedEdge& child_inner_to_outer_edge,
        const CompressedEdge& child_compressed_edge);
    AltTreePruneResult prune_upward_path_stopping_before(AltTreeNode* prune_parent, bool back);
    AltTreePruneResult prune_upward_back_edge_path_stopping_before(AltTreeNode* prune_parent);
    const AltTreeNode* find_root() const;
    bool tree_equal(const pm::AltTreeNode& other) const;
    std::vector<pm::AltTreeNode*> all_nodes_in_tree();
};

}  // namespace pm

#endif  // PYMATCHING2_ALTERNATING_TREE_H
