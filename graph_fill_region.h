#ifndef PYMATCHING2_GRAPH_FILL_REGION_H
#define PYMATCHING2_GRAPH_FILL_REGION_H

#include "alternating_tree.h"
#include "region_edge.h"


namespace pm {

    class AltTreeNode;

    struct Match {
        pm::GraphFillRegion* region;
        pm::CompressedEdge edge;
        Match(pm::GraphFillRegion* region, pm::CompressedEdge edge);
        Match();
    };

    class GraphFillRegion {
    public:
        // Maybe one field for blossom_parent or alt_tree_node eventually by using union?
        GraphFillRegion* blossom_parent;
        pm::AltTreeNode* alt_tree_node;
        pm::Varying32 radius;

        pm::Match match;

        std::vector<pm::RegionEdge> blossom_children;
        std::vector<pm::DetectorNode*> shell_area;

        GraphFillRegion();
        bool tree_equal(const pm::GraphFillRegion& other) const;

        const pm::GraphFillRegion* top_region() const;
        void add_match(pm::GraphFillRegion* match, const pm::CompressedEdge& edge);

        bool operator==(const GraphFillRegion &rhs) const;

        bool operator!=(const GraphFillRegion &rhs) const;
    };

}


#endif //PYMATCHING2_GRAPH_FILL_REGION_H
