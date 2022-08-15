#ifndef PYMATCHING2_GRAPH_FLOODER_H
#define PYMATCHING2_GRAPH_FLOODER_H

#include <queue>

#include "pymatching/events.h"
#include "pymatching/graph.h"
#include "pymatching/region_edge.h"
#include "pymatching/bit_bucket_queue.h"

namespace pm {

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
    bit_bucket_queue<false> queue;

    /// This counter is used to give each tentative event a unique index. The index is also written
    /// to objects affected by the event. For an event to be valid, its vid must match the marked
    /// vids on the objects it affects. This allows you to invalidate all events affecting an
    /// object by incrementing its vids.
    uint64_t next_event_vid;

    GraphFlooder(MatchingGraph graph);
    GraphFlooder(GraphFlooder&&) noexcept;
    void create_region(DetectorNode* node);
    MwpmEvent next_event();
    void set_region_growing(pm::GraphFillRegion& region);
    void set_region_frozen(pm::GraphFillRegion& region);
    void set_region_shrinking(pm::GraphFillRegion& region);
    GraphFillRegion* create_blossom(std::vector<RegionEdge>& contained_regions);
    void schedule_tentative_shrink_event(GraphFillRegion& region);
    void reschedule_events_at_detector_node(DetectorNode& detector_node);
    void do_region_created_at_empty_detector_node(GraphFillRegion& region, DetectorNode& detector_node);
    void do_region_arriving_at_empty_detector_node(
        GraphFillRegion& region, DetectorNode& empty_node, DetectorNode& from_node, size_t neighbor_index);
    MwpmEvent do_region_shrinking(const TentativeEventData_LookAtShrinkingRegion& event);
    pm::MwpmEvent do_neighbor_interaction(
        DetectorNode &src,
        size_t src_to_dst_index,
        DetectorNode &dst,
        size_t dst_to_src_index);
    pm::MwpmEvent do_region_hit_boundary_interaction(DetectorNode &node);
    static MwpmEvent do_degenerate_implosion(const GraphFillRegion& region);
    static MwpmEvent do_blossom_shattering(GraphFillRegion& region);
    bool dequeue_decision(pm::TentativeEvent ev);
    std::pair<size_t, pm::cumulative_time_int> find_next_event_at_node_returning_neighbor_index_and_time(DetectorNode &detector_node) const;
    pm::MwpmEvent do_look_at_node_event(DetectorNode &node);

    pm::TentativeEvent dequeue_valid();
    pm::MwpmEvent do_valid_tentative_event_returning_mwpm_event(TentativeEvent tentative_event);
};

}  // namespace pm

#endif  // PYMATCHING2_GRAPH_FLOODER_H
