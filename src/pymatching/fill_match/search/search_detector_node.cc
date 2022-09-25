
#include "search_detector_node.h"

size_t pm::SearchDetectorNode::index_of_neighbor(SearchDetectorNode *target) const {
    for (size_t k = 0; k < neighbors.size(); k++) {
        if (neighbors[k] == target) {
            return k;
        }
    }
    throw std::invalid_argument("Failed to find neighbor.");
}

void pm::SearchDetectorNode::reset() {
    reached_from_source = nullptr;
    index_of_predecessor = SIZE_MAX;
    node_event_tracker.clear();
}

