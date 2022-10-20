#ifndef PYMATCHING2_SEARCH_GRAPH_H
#define PYMATCHING2_SEARCH_GRAPH_H

#include "pymatching/sparse_blossom/search/search_detector_node.h"
#include "pymatching/sparse_blossom/tracker/queued_event_tracker.h"

namespace pm {

class SearchGraph {
   public:
    std::vector<SearchDetectorNode> nodes;
    size_t num_nodes;

    SearchGraph();
    explicit SearchGraph(size_t num_nodes);
    SearchGraph(SearchGraph&& graph) noexcept;
    void add_edge(size_t u, size_t v, signed_weight_int weight, const std::vector<size_t>& observables);
    void add_boundary_edge(size_t u, signed_weight_int weight, const std::vector<size_t>& observables);
};
}  // namespace pm

#endif  // PYMATCHING2_SEARCH_GRAPH_H
