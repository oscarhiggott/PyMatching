#ifndef PYMATCHING2_SEARCH_FLOODER_H
#define PYMATCHING2_SEARCH_FLOODER_H

#include "pymatching/fill_match/search/search_graph.h"
#include "pymatching/fill_match/tracker/radix_heap_queue.h"

namespace pm {

enum TargetType : uint8_t { DETECTOR_NODE, BOUNDARY, NO_TARGET };

/// The edge on which two search regions collided
struct CollisionEdge {
    SearchDetectorNode* collision_node;
    size_t neighbor_index;
};

class SearchFlooder {
   public:
    SearchFlooder();
    explicit SearchFlooder(SearchGraph graph);
    SearchFlooder(SearchFlooder&& other) noexcept;
    SearchGraph graph;
    pm::radix_heap_queue<false> queue;
    /// The reached_nodes are the nodes that need to be reset after each search completes
    std::vector<SearchDetectorNode*> reached_nodes;
    /// The type of target for the search from a detection event, either another detection event or the boundary.
    TargetType target_type;
    void reschedule_events_at_search_detector_node(SearchDetectorNode& detector_node);
    std::pair<size_t, cumulative_time_int> find_next_event_at_node_returning_neighbor_index_and_time(
        const SearchDetectorNode& detector_node) const;
    void do_search_starting_at_empty_search_detector_node(SearchDetectorNode* src);
    void do_search_exploring_empty_detector_node(SearchDetectorNode& empty_node, size_t empty_to_from_index);
    CollisionEdge do_look_at_node_event(SearchDetectorNode& node);
    CollisionEdge run_until_collision(SearchDetectorNode* src, SearchDetectorNode* dst);
    void trace_back_path_from_node(
        SearchDetectorNode* detector_node, uint8_t* obs_begin_ptr, pm::total_weight_int& weight);
    void trace_back_path_from_collision_edge(
        CollisionEdge collision_edge, uint8_t* obs_begin_ptr, pm::total_weight_int& weight);
    void reset_graph();
    void reset();
};

}  // namespace pm

#endif  // PYMATCHING2_SEARCH_FLOODER_H
