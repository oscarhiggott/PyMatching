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
        pm::TentativeEvent* shrink_event;
        pm::Match match;

        std::vector<pm::RegionEdge> blossom_children;
        std::vector<pm::DetectorNode*> shell_area;



        GraphFillRegion();
        bool tree_equal(const pm::GraphFillRegion& other) const;

        GraphFillRegion * top_region() const;
        void add_match(pm::GraphFillRegion* match, const pm::CompressedEdge& edge);

        template<typename Callable>
        void do_op_for_each_node_in_total_area(const Callable& func);

        bool operator==(const GraphFillRegion &rhs) const;

        bool operator!=(const GraphFillRegion &rhs) const;
    };

    template<typename Callable>
    inline void pm::GraphFillRegion::do_op_for_each_node_in_total_area(const Callable& func) {
        for (size_t i = 0; i < shell_area.size(); i++){
            func(shell_area[shell_area.size() - i - 1]);
        }
        for (auto& child: blossom_children) {
            child.region->do_op_for_each_node_in_total_area(func);
        }
    }

}


#endif //PYMATCHING2_GRAPH_FILL_REGION_H
