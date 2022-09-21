#ifndef PYMATCHING2_SEARCH_GRAPH_H
#define PYMATCHING2_SEARCH_GRAPH_H

#include "pymatching/fill_match/tracker/queued_event_tracker.h"

namespace pm {

    class SearchDetectorNode {
    public:
        SearchDetectorNode()
            : reached_from_source(nullptr),
              index_of_predecessor(SIZE_MAX) {
        }

        /// The SearchDetectorNode that this node was reached from in the Dijkstra search
        SearchDetectorNode* reached_from_source;

        /// `index_of_predecessor` is the index in `neighbors` of the neighboring detector node that this node was
        /// reached from in the Dijkstra search.
        size_t index_of_predecessor;

        /// Manages the next "look at me" event for the node
        QueuedEventTracker node_event_tracker;

        /// == Permanent fields used to define the structure of the graph. ==
        std::vector<DetectorNode*> neighbors;       /// The node's neighbors.
        std::vector<weight_int> neighbor_weights;   /// Distance crossed by the edge to each neighbor.
        std::vector<std::vector<size_t>> neighbor_observable_indices;  /// Indices of observables crossed by the edge to each neighbor.
    };
}


#endif //PYMATCHING2_SEARCH_GRAPH_H
