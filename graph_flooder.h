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
        inline bool operator() (TentativeEvent a, TentativeEvent b) {
            return a > b;
        }
    };

    class GraphFlooder {
    public:
        /// The graph of detector nodes that is being flooded.
        MatchingGraph graph;
        /// Tracks the next thing that will occur as flooding proceeds.
        /// The events are "tentative" because processing an event may remove another,
        /// for example if a region stops growing due to colliding with another region
        /// then the growing region will no longer reach other nodes even if those
        /// events were scheduled in this queue before the region collision was processed.
        ///
        /// Events are ordered by time; by when they will occur in a timeline.
        std::priority_queue<TentativeEvent, std::vector<TentativeEvent>, Compare> raw_queue;
        /// Tracks where the flooder is in the timeline of events to process.
        ///
        /// As time advances, regions grow and shrink and collide and change.
        /// Because regions grow at a finite rate, collisions and other events
        /// occur at specific times. This allows them to be processed in the
        /// correct order.
        time_int time;
        size_t next_validation_index;

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

//        // TODO: Put in class called bucket queue. Template with number of buckets. Switch to heapq if number of buckets is big.
//        // Use standard library heap
//        std::vector<TentativeEvent*> bucket_queue[MAX_TIME];
    };

}

#endif //PYMATCHING2_GRAPH_FLOODER_H
