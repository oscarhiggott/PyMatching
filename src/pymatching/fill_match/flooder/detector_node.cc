#include "pymatching/fill_match/flooder/detector_node.h"

#include "pymatching/fill_match/flooder/graph_fill_region.h"

namespace pm {

Varying32 DetectorNode::total_radius() const {
    Varying32 result{0};
    if (reached_from_source != nullptr) {
        auto r = reached_from_source->region_that_arrived;
        while (r != nullptr) {
            result.inplace_freeze_then_add(r->radius);
            r = r->blossom_parent;
        }
    }
    return result;
}

Varying32 DetectorNode::local_radius() const {
    return total_radius() - distance_from_source;
}

bool DetectorNode::has_same_owner_as(const DetectorNode &other) const {
    if ((region_that_arrived == nullptr) != (other.region_that_arrived == nullptr))
        return false;
    if (region_that_arrived == other.region_that_arrived)
        return true;
    return top_region() == other.top_region();
}

GraphFillRegion *DetectorNode::top_region() const {
    if (!region_that_arrived)
        return nullptr;
    return region_that_arrived->top_region();
}

void DetectorNode::reset() {
    observables_crossed_from_source = 0;
    reached_from_source = nullptr;
    distance_from_source = 0;
    region_that_arrived = nullptr;
    node_event_tracker.clear();
}

size_t DetectorNode::index_of_neighbor(DetectorNode *target) const {
    for (size_t k = 0; k < neighbors.size(); k++) {
        if (neighbors[k] == target) {
            return k;
        }
    }
    throw std::invalid_argument("Failed to find neighbor.");
}

}  // namespace pm
