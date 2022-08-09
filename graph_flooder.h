#ifndef PYMATCHING2_GRAPH_FLOODER_H
#define PYMATCHING2_GRAPH_FLOODER_H

#include "graph.h"
#include "events.h"
#include "region_edge.h"
#include<queue>


namespace pm{

//    constexpr int MAX_TIME = 100000; // twice maximum edge weight

    class Compare{
    public:
        bool operator() (TentativeEvent* a, TentativeEvent* b) {
            return *a > *b;
        }
    };

    class GraphFlooder {
    public:
        MatchingGraph graph;
        std::priority_queue<TentativeEvent*, std::vector<TentativeEvent*>, Compare> queue;
        time_int time;
        explicit GraphFlooder(MatchingGraph& graph);
        GraphFlooder(GraphFlooder&&) noexcept;
        void create_region(DetectorNode* node);
        MwpmEvent next_event();
        void set_region_growing(pm::GraphFillRegion& region);
        void set_region_frozen(pm::GraphFillRegion& region);
        void set_region_shrinking(pm::GraphFillRegion& region);
        GraphFillRegion* create_blossom(std::vector<RegionEdge>& contained_regions);
        void schedule_tentative_neighbor_interaction_event(
                DetectorNode* detector_node_1,
                size_t schedule_list_index_1,
                DetectorNode* detector_node_2,
                size_t schedule_list_index_2,
                time_int event_time
        );
        void schedule_tentative_shrink_event(GraphFillRegion& region);
        void reschedule_events_for_region(GraphFillRegion& region);
        void reschedule_events_at_detector_node(DetectorNode& detector_node);
        void do_region_created_at_empty_detector_node(GraphFillRegion& region, DetectorNode& detector_node);
        void do_region_arriving_at_empty_detector_node(GraphFillRegion& region, DetectorNode& empty_node,
                                                       DetectorNode& from_node, size_t neighbor_index);
        MwpmEvent do_region_shrinking(const TentativeRegionShrinkEventData& event);
        MwpmEvent do_neighbor_interaction(const TentativeNeighborInteractionEventData& event);
        static MwpmEvent do_region_hit_boundary_interaction(const TentativeNeighborInteractionEventData& event);
        static MwpmEvent do_degenerate_implosion(const GraphFillRegion& region);
        static MwpmEvent do_blossom_shattering(GraphFillRegion& region);

//        // Put in class called bucket queue. Template with number of buckets. Switch to heapq if number of buckets is big.
//        // Use standard library heap
//        std::vector<TentativeEvent*> bucket_queue[MAX_TIME];
    };

}

#endif //PYMATCHING2_GRAPH_FLOODER_H
