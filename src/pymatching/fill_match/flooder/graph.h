#ifndef PYMATCHING2_GRAPH_H
#define PYMATCHING2_GRAPH_H

#include <vector>

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
    size_t num_nodes;

    MatchingGraph();
    explicit MatchingGraph(size_t num_nodes);
    MatchingGraph(MatchingGraph&& graph) noexcept;
    void add_edge(size_t u, size_t v, weight_int weight, obs_int observables);
    void add_edge(size_t u, size_t v, weight_int weight, const std::vector<size_t>& observables);
    void add_boundary_edge(size_t u, weight_int weight, obs_int observables);
    void add_boundary_edge(size_t u, weight_int weight, const std::vector<size_t>& observables);
};

}  // namespace pm

#endif  // PYMATCHING2_GRAPH_H
