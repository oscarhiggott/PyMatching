#ifndef PYMATCHING2_ALTERNATING_TREE_H
#define PYMATCHING2_ALTERNATING_TREE_H


#include "compressed_edge.h"
#include "region_path.h"
#include "graph_fill_region.h"

namespace pm{

    class AltTreeNode;
    class GraphFillRegion;


    struct AltTreeEdge{
        AltTreeNode* alt_tree_node;
        CompressedEdge edge;
        AltTreeEdge();
        AltTreeEdge(AltTreeNode* alt_tree_node, const CompressedEdge& edge);
    };

// Equality and string


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
        std::vector<GraphFillRegion> shatter_into_matches();
        AltTreeNode most_recent_common_ancestor(const AltTreeNode& other);
        bool in_same_tree_as(const AltTreeNode& other);
        void add_child(const AltTreeEdge& child);
        AltTreeNode* make_child(GraphFillRegion* inner_region, GraphFillRegion* outer_region,
                                const CompressedEdge& inner_outer_edge, const CompressedEdge& child_edge);
        AltTreePruneResult prune_upward_stopping_before(AltTreeNode* prune_parent);
        AltTreeNode* find_root();
    };

}

#endif //PYMATCHING2_ALTERNATING_TREE_H
