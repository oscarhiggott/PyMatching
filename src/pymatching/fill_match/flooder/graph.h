#ifndef PYMATCHING2_GRAPH_H
#define PYMATCHING2_GRAPH_H

#include <vector>
#include <set>

#include "pymatching/fill_match/flooder/detector_node.h"
#include "pymatching/fill_match/flooder_matcher_interop/varying.h"
#include "pymatching/fill_match/tracker/flood_check_event.h"
#include "pymatching/fill_match/tracker/queued_event_tracker.h"

namespace pm {

class GraphFillRegion;

/// A collection of detector nodes. It's expected that all detector nodes in the graph
/// will only refer to other detector nodes within the same graph.
class MatchingGraph {
   public:
    std::vector<DetectorNode> nodes;
    /// These are the detection events that would occur if an error occurred on every edge with a negative weight
    std::set<size_t> negative_weight_detection_events_set;
    /// These are the observables that would be flipped if an error occurred on every edge with a negative weight
    std::set<size_t> negative_weight_observables_set;
    /// The sum of the negative edge weights. This number is negative, rather than the absolute value.
    pm::total_weight_int negative_weight_sum;
    size_t num_nodes;
    size_t num_observables;

    MatchingGraph();
    explicit MatchingGraph(size_t num_nodes, size_t num_observables);
    MatchingGraph(MatchingGraph&& graph) noexcept;
    void add_edge(size_t u, size_t v, signed_weight_int weight, const std::vector<size_t>& observables);
    void add_boundary_edge(size_t u, signed_weight_int weight, const std::vector<size_t>& observables);
    void update_negative_weight_observables(const std::vector<size_t>& observables);
    void update_negative_weight_detection_events(size_t node_id);
};

}  // namespace pm

#endif  // PYMATCHING2_GRAPH_H
