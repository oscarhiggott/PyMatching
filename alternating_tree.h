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
    };

// Equality and string


    struct AltTreePruneResult {
        std::vector<AltTreeEdge> orphans;
        std::vector<RegionEdge> pruned_path_regions;
    };

    class AltTreeNode{
    public:
        GraphFillRegion* inner_region;
        GraphFillRegion* outer_region;
        AltTreeEdge parent;
        std::vector<AltTreeEdge> children; // Maybe make linked list?
        CompressedEdge inner_to_outer_edge;

        std::vector<GraphFillRegion> shatter_into_matches();
        AltTreeNode most_recent_common_ancestor(AltTreeNode other);
        bool in_same_tree_as(AltTreeNode other);
        void add_child(AltTreeEdge child);
        AltTreeNode* make_child(GraphFillRegion* inner_region, GraphFillRegion* outer_region,
                                CompressedEdge inner_outer_edge, CompressedEdge child_edge);
        AltTreePruneResult prune_upward_stopping_before(AltTreeNode* prune_parent);
    };

}

#endif //PYMATCHING2_ALTERNATING_TREE_H
