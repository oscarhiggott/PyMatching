#ifndef PYMATCHING2_GRAPH_FLOODER_H
#define PYMATCHING2_GRAPH_FLOODER_H

#include "graph.h"
#include "events.h"
#include "region_edge.h"

// add namespace
//
//constexpr int MAX_TIME = 100000; // twice maximum edge weight
//
//
//class GraphFlooder {
//public:
//    Graph graph;
//    pm::time_int time;
//    MwpmEvent next_event();
//    void set_region_growth(GraphFillRegion& region);
//    GraphFillRegion* create_blossom(const RegionPath& contained_regions);
//    void schedule_tentative_neighbor_interaction_event(
//            DetectorNode& detector_node_1,
//            int schedule_list_index_1,
//            DetectorNode& detector_node_2,
//            int schedule_list_index_2
//    );
//    void schedule_tentative_shrink_event(GraphFillRegion& region);
//    void reschedule_events_for_region(GraphFillRegion& region);
//    void reschedule_events_at_detector_node(DetectorNode& detector_node);
//    void do_region_arriving_at_empty_detector_node(const GraphFillRegion& region, DetectorNode& detector_node);
//    MwpmEvent do_region_shrinking(const TentativeRegionShrinkEventData& event); // Use std::optional?
//    MwpmEvent do_neighbor_interaction(const TentativeNeighborInteractionEvent& event); // Use std::optional?
//    MwpmEvent do_region_hit_boundary_interaction(const TentativeNeighborInteractionEvent& event);
//    MwpmEvent do_degenerate_implosion(const GraphFillRegion& region);
//    MwpmEvent do_blossom_implosioin(GraphFillRegion& region);
//
//    // Delete these:
//    int next_event_id;
//    int next_region_id;
//    // Put in class called bucket queue. Template with number of buckets. Switch to heapq if number of buckets is big.
//    // Use standard library heap
//    std::vector<TentativeEvent*> bucket_queue[MAX_TIME];
//};

#endif //PYMATCHING2_GRAPH_FLOODER_H
