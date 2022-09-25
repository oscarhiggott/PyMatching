
#include "search_detector_node.h"

size_t pm::SearchDetectorNode::index_of_neighbor(SearchDetectorNode *target) const {
    for (size_t k = 0; k < neighbors.size(); k++) {
        if (neighbors[k] == target) {
            return k;
        }
    }
    throw std::invalid_argument("Failed to find neighbor.");
}

