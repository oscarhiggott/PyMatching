#ifndef PYMATCHING2_GRAPH_FILL_REGION_H
#define PYMATCHING2_GRAPH_FILL_REGION_H

#include "graph.h"
#include "alternating_tree.h"

namespace pm {
    class GraphFillRegion {
    public:
        // Maybe one field for blossom_parent or alt_tree_node eventually by using union?
        GraphFillRegion* blossom_parent;
        pm::AltTreeNode* alt_tree_node;
        pm::Varying32 radius;

        std::vector<GraphFillRegion*> blossom_children;
        std::vector<pm::DetectorNode*> shell_area;
    };
}


#endif //PYMATCHING2_GRAPH_FILL_REGION_H
