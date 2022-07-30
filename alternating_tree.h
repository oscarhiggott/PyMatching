#ifndef PYMATCHING2_ALTERNATING_TREE_H
#define PYMATCHING2_ALTERNATING_TREE_H


#include "compressed_edge.h"
#include "region_edge.h"
#include "graph_fill_region.h"

namespace pm{

    class AltTreeNode;
    class GraphFillRegion;


    struct AltTreeEdge{
        AltTreeNode* alt_tree_node;
        pm::CompressedEdge edge;
        AltTreeEdge();

        bool operator==(const AltTreeEdge &rhs) const;

        bool operator!=(const AltTreeEdge &rhs) const;

        AltTreeEdge(AltTreeNode* alt_tree_node, const CompressedEdge& edge);
    };

    template<class T, class UnaryPredicate>
    bool unstable_erase(std::vector<T>& vec, UnaryPredicate pred);

    template<class T, class UnaryPredicate>
    inline bool unstable_erase(std::vector<T>& vec, UnaryPredicate pred) {
        auto res = std::find_if(vec.begin(), vec.end(), pred);
        if (res == vec.end())
            return false;
        if (vec.size() > 1){
            *res = std::move(vec.back());
            vec.pop_back();
            return true;
        }
        vec.clear();
        return true;
    }

    struct AltTreePruneResult {
        std::vector<AltTreeEdge> orphans;
        std::vector<RegionEdge> pruned_path_regions;
    };

    class AltTreeNode {
    public:
        GraphFillRegion* inner_region;
        GraphFillRegion* outer_region;
        CompressedEdge inner_to_outer_edge;
        AltTreeEdge parent;
        std::vector<AltTreeEdge> children; // Maybe make linked list?
        AltTreeNode();
        AltTreeNode(GraphFillRegion* inner_region, GraphFillRegion* outer_region,
                    const CompressedEdge& inner_to_outer_edge, const AltTreeEdge& parent,
                    std::vector<AltTreeEdge>  children);
        AltTreeNode(pm::GraphFillRegion* inner_region, pm::GraphFillRegion* outer_region,
                    const pm::CompressedEdge& inner_to_outer_edge);
        AltTreeNode(pm::GraphFillRegion* outer_region);
        ~AltTreeNode();

        bool operator==(const AltTreeNode &rhs) const;

        bool operator!=(const AltTreeNode &rhs) const;

        void become_root();
        std::vector<GraphFillRegion> shatter_into_matches();
        AltTreeNode most_recent_common_ancestor(const AltTreeNode& other);
        bool in_same_tree_as(const AltTreeNode& other);
        void add_child(const AltTreeEdge& child);
        AltTreeNode* make_child(GraphFillRegion* child_inner_region, GraphFillRegion* child_outer_region,
                                const CompressedEdge& child_inner_to_outer_edge, const CompressedEdge& child_compressed_edge);
        AltTreePruneResult prune_upward_stopping_before(AltTreeNode* prune_parent);
        const AltTreeNode* find_root() const;
        bool tree_equal(const pm::AltTreeNode& other) const;
        std::vector<pm::AltTreeNode*> all_nodes_in_tree();
    };

}

#endif //PYMATCHING2_ALTERNATING_TREE_H
