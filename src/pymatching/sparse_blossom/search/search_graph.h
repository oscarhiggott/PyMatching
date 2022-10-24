#ifndef PYMATCHING2_SEARCH_GRAPH_H
#define PYMATCHING2_SEARCH_GRAPH_H

#include "pymatching/sparse_blossom/search/search_detector_node.h"
#include "pymatching/sparse_blossom/tracker/queued_event_tracker.h"

namespace pm {

/// The edge on which two search regions collided
struct SearchGraphEdge {
    SearchDetectorNode* detector_node;
    size_t neighbor_index;
};

class SearchGraph {
   public:
    std::vector<SearchDetectorNode> nodes;
    size_t num_nodes;
    std::vector<SearchGraphEdge> negative_weight_edges;

    SearchGraph();
    explicit SearchGraph(size_t num_nodes);
    SearchGraph(SearchGraph&& graph) noexcept;
    void add_edge(size_t u, size_t v, signed_weight_int weight, const std::vector<size_t>& observables);
    void add_boundary_edge(size_t u, signed_weight_int weight, const std::vector<size_t>& observables);
};
}  // namespace pm

#endif  // PYMATCHING2_SEARCH_GRAPH_H
