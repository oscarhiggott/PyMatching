#include "pymatching/fill_match/flooder/detector_node.h"

#include "pymatching/fill_match/flooder/graph_fill_region.h"

namespace pm {

int32_t DetectorNode::compute_wrapped_radius() const {
    if (reached_from_source == nullptr) {
        return 0;
    }
    int32_t total = 0;
    auto r = region_that_arrived;
    while (r != region_that_arrived_top) {
        total += r->radius.y_intercept();
        r = r->blossom_parent;
    }
    return total - radius_of_arrival;
}

cumulative_time_int DetectorNode::compute_wrapped_radius_within_layer_at_time(GraphFillRegion *target_layer, cumulative_time_int t) const {
    int32_t radius_past_region_that_arrived = 0;
    auto r = region_that_arrived;
    while (r != nullptr) {
        radius_past_region_that_arrived += r->radius.get_distance_at_time(t);
        if (r == target_layer) {
            break;
        }
        r = r->blossom_parent;
    }
    return radius_past_region_that_arrived - radius_of_arrival;
}

void DetectorNode::reset() {
    observables_crossed_from_source = 0;
    reached_from_source = nullptr;
    radius_of_arrival = 0;
    region_that_arrived = nullptr;
    region_that_arrived_top = nullptr;
    wrapped_radius_cached = 0;
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

GraphFillRegion *DetectorNode::heir_region_on_shatter() const {
    GraphFillRegion *r = region_that_arrived;
    while (true) {
        GraphFillRegion *p = r->blossom_parent;
        if (p == region_that_arrived_top) {
            return r;
        }
        r = p;
    }
}

}  // namespace pm
