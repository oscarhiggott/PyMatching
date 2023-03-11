// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "pymatching/sparse_blossom/flooder/graph_flooder.h"

#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"
#include "pymatching/sparse_blossom/matcher/alternating_tree.h"

using namespace pm;

GraphFlooder::GraphFlooder(MatchingGraph graph)
    : graph(std::move(graph)), negative_weight_obs_mask(0), negative_weight_sum(0) {
}

GraphFlooder::GraphFlooder(GraphFlooder &&flooder) noexcept
    : graph(std::move(flooder.graph)),
      queue(std::move(flooder.queue)),
      region_arena(std::move(flooder.region_arena)),
      match_edges(std::move(flooder.match_edges)),
      negative_weight_detection_events(std::move(flooder.negative_weight_detection_events)),
      negative_weight_observables(std::move(flooder.negative_weight_observables)),
      negative_weight_obs_mask(flooder.negative_weight_obs_mask),
      negative_weight_sum(flooder.negative_weight_sum) {
}

void GraphFlooder::do_region_created_at_empty_detector_node(GraphFillRegion &region, DetectorNode &detector_node) {
    detector_node.reached_from_source = &detector_node;
    detector_node.radius_of_arrival = 0;
    detector_node.region_that_arrived = &region;
    detector_node.region_that_arrived_top = &region;
    detector_node.wrapped_radius_cached = 0;
    region.shell_area.push_back(&detector_node);
    reschedule_events_at_detector_node(detector_node);
}

std::pair<size_t, cumulative_time_int> find_next_event_at_node_not_occupied_by_growing_top_region(
    const DetectorNode &detector_node, VaryingCT rad1) {
    cumulative_time_int best_time = std::numeric_limits<cumulative_time_int>::max();
    size_t best_neighbor = SIZE_MAX;

    size_t start = 0;
    if (!detector_node.neighbors.empty() && detector_node.neighbors[0] == nullptr)
        start++;

    // Handle non-boundary neighbors.
    for (size_t i = start; i < detector_node.neighbors.size(); i++) {
        auto weight = detector_node.neighbor_weights[i];

        auto neighbor = detector_node.neighbors[i];

        auto rad2 = neighbor->local_radius();

        if (rad2.is_growing()) {
            auto collision_time = weight - rad1.y_intercept() - rad2.y_intercept();
            if (collision_time < best_time) {
                best_time = collision_time;
                best_neighbor = i;
            }
        }
    }
    return {best_neighbor, best_time};
}

std::pair<size_t, cumulative_time_int> find_next_event_at_node_occupied_by_growing_top_region(
    const DetectorNode &detector_node, const VaryingCT &rad1) {
    cumulative_time_int best_time = std::numeric_limits<cumulative_time_int>::max();
    size_t best_neighbor = SIZE_MAX;
    size_t start = 0;
    if (!detector_node.neighbors.empty() && detector_node.neighbors[0] == nullptr) {
        // Growing towards boundary
        auto weight = detector_node.neighbor_weights[0];
        auto collision_time = weight - rad1.y_intercept();
        if (collision_time < best_time) {
            best_time = collision_time;
            best_neighbor = 0;
        }
        start++;
    }

    // Handle non-boundary neighbors.
    for (size_t i = start; i < detector_node.neighbors.size(); i++) {
        auto weight = detector_node.neighbor_weights[i];

        auto neighbor = detector_node.neighbors[i];
        if (detector_node.has_same_owner_as(*neighbor)) {
            continue;
        }
        auto rad2 = neighbor->local_radius();
        if (rad2.is_shrinking()) {
            continue;
        }

        auto collision_time = weight - rad1.y_intercept() - rad2.y_intercept();
        if (rad2.is_growing()) {
            collision_time >>= 1;
        }
        if (collision_time < best_time) {
            best_time = collision_time;
            best_neighbor = i;
        }
    }
    return {best_neighbor, best_time};
}

std::pair<size_t, cumulative_time_int> GraphFlooder::find_next_event_at_node_returning_neighbor_index_and_time(
    const DetectorNode &detector_node) const {
    auto rad1 = detector_node.local_radius();

    if (rad1.is_growing()) {
        return find_next_event_at_node_occupied_by_growing_top_region(detector_node, rad1);
    } else {
        return find_next_event_at_node_not_occupied_by_growing_top_region(detector_node, rad1);
    }
}

void GraphFlooder::reschedule_events_at_detector_node(DetectorNode &detector_node) {
    auto x = find_next_event_at_node_returning_neighbor_index_and_time(detector_node);
    if (x.first == SIZE_MAX) {
        detector_node.node_event_tracker.set_no_desired_event();
    } else {
        detector_node.node_event_tracker.set_desired_event(
            {
                &detector_node,
                cyclic_time_int{x.second},
            },
            queue);
    }
}

void GraphFlooder::schedule_tentative_shrink_event(GraphFillRegion &region) {
    cumulative_time_int t;
    if (region.shell_area.empty()) {
        t = region.radius.time_of_x_intercept_for_shrinking();
    } else {
        t = region.shell_area.back()->local_radius().time_of_x_intercept_for_shrinking();
    }
    region.shrink_event_tracker.set_desired_event(
        {
            &region,
            cyclic_time_int{t},
        },
        queue);
}

void GraphFlooder::do_region_arriving_at_empty_detector_node(
    GraphFillRegion &region, DetectorNode &empty_node, const DetectorNode &from_node, size_t from_to_empty_index) {
    empty_node.observables_crossed_from_source =
        (from_node.observables_crossed_from_source ^ from_node.neighbor_observables[from_to_empty_index]);
    empty_node.reached_from_source = from_node.reached_from_source;
    empty_node.radius_of_arrival = region.radius.get_distance_at_time(queue.cur_time);
    empty_node.region_that_arrived = &region;
    empty_node.region_that_arrived_top = region.blossom_parent_top;
    empty_node.wrapped_radius_cached = empty_node.compute_wrapped_radius();
    region.shell_area.push_back(&empty_node);
    reschedule_events_at_detector_node(empty_node);
}

MwpmEvent GraphFlooder::do_region_shrinking(GraphFillRegion &region) {
    if (region.shell_area.empty()) {
        return do_blossom_shattering(region);
    } else if (region.shell_area.size() == 1 && region.blossom_children.empty()) {
        return do_degenerate_implosion(region);
    } else {
        auto leaving_node = region.shell_area.back();
        region.shell_area.pop_back();
        leaving_node->region_that_arrived = nullptr;
        leaving_node->region_that_arrived_top = nullptr;
        leaving_node->wrapped_radius_cached = 0;
        leaving_node->reached_from_source = nullptr;
        leaving_node->radius_of_arrival = 0;
        leaving_node->observables_crossed_from_source = 0;
        reschedule_events_at_detector_node(*leaving_node);
        schedule_tentative_shrink_event(region);
        return MwpmEvent::no_event();
    }
}

MwpmEvent GraphFlooder::do_neighbor_interaction(DetectorNode &src, size_t src_to_dst_index, DetectorNode &dst) {
    // First check if one region is moving into an empty location
    if (src.region_that_arrived && !dst.region_that_arrived) {
        do_region_arriving_at_empty_detector_node(*src.region_that_arrived_top, dst, src, src_to_dst_index);
        return MwpmEvent::no_event();
    } else if (dst.region_that_arrived && !src.region_that_arrived) {
        do_region_arriving_at_empty_detector_node(*dst.region_that_arrived_top, src, dst, dst.index_of_neighbor(&src));
        return MwpmEvent::no_event();
    } else {
        // Two regions colliding
        return RegionHitRegionEventData{
            src.region_that_arrived_top,
            dst.region_that_arrived_top,
            CompressedEdge{
                src.reached_from_source,
                dst.reached_from_source,
                src.observables_crossed_from_source ^ dst.observables_crossed_from_source ^
                    src.neighbor_observables[src_to_dst_index]},
        };
    }
}

MwpmEvent GraphFlooder::do_region_hit_boundary_interaction(DetectorNode &node) {
    return RegionHitBoundaryEventData{
        node.region_that_arrived_top,
        CompressedEdge{
            node.reached_from_source, nullptr, node.observables_crossed_from_source ^ node.neighbor_observables[0]}};
}

MwpmEvent GraphFlooder::do_degenerate_implosion(const GraphFillRegion &region) {
    return RegionHitRegionEventData{
        region.alt_tree_node->parent.alt_tree_node->outer_region,
        region.alt_tree_node->outer_region,
        CompressedEdge{
            region.alt_tree_node->parent.edge.loc_to,
            region.alt_tree_node->inner_to_outer_edge.loc_to,
            region.alt_tree_node->inner_to_outer_edge.obs_mask ^ region.alt_tree_node->parent.edge.obs_mask}};
}

MwpmEvent GraphFlooder::do_blossom_shattering(GraphFillRegion &region) {
    return BlossomShatterEventData{
        &region,
        region.alt_tree_node->parent.edge.loc_from->heir_region_on_shatter(),
        region.alt_tree_node->inner_to_outer_edge.loc_from->heir_region_on_shatter()};
}

GraphFillRegion *GraphFlooder::create_blossom(std::vector<RegionEdge> &contained_regions) {
    auto blossom_region = region_arena.alloc_default_constructed();
    blossom_region->radius = VaryingCT::growing_varying_with_zero_distance_at_time(queue.cur_time);
    blossom_region->blossom_children = std::move(contained_regions);
    for (auto &region_edge : blossom_region->blossom_children) {
        region_edge.region->radius = region_edge.region->radius.then_frozen_at_time(queue.cur_time);
        region_edge.region->wrap_into_blossom(blossom_region);
        region_edge.region->shrink_event_tracker.set_no_desired_event();
    }

    blossom_region->do_op_for_each_node_in_total_area([this](DetectorNode *n) {
        reschedule_events_at_detector_node(*n);
    });

    return blossom_region;
}

bool GraphFlooder::dequeue_decision(FloodCheckEvent ev) {
    switch (ev.tentative_event_type) {
        case FloodCheckEventType::LOOK_AT_NODE: {
            auto &node = *ev.data_look_at_node;
            return node.node_event_tracker.dequeue_decision(ev, queue);
        }
        case FloodCheckEventType::LOOK_AT_SHRINKING_REGION: {
            auto &region = *ev.data_look_at_shrinking_region;
            return region.shrink_event_tracker.dequeue_decision(ev, queue);
        }
        case FloodCheckEventType::NO_FLOOD_CHECK_EVENT:
            return true;
        default:
            throw std::invalid_argument("Unrecognized event type.");
    }
}

FloodCheckEvent GraphFlooder::dequeue_valid() {
    while (true) {
        FloodCheckEvent ev = queue.dequeue();
        if (dequeue_decision(ev)) {
            return ev;
        }
    }
}

void GraphFlooder::set_region_growing(GraphFillRegion &region) {
    region.radius = region.radius.then_growing_at_time(queue.cur_time);

    // No shrinking event can occur while growing.
    region.shrink_event_tracker.set_no_desired_event();

    // Node events can occur while growing, and events in the queue may occur sooner than
    // previously scheduled. Therefore, we must reschedule all the nodes.
    region.do_op_for_each_node_in_total_area([this](DetectorNode *n) {
        reschedule_events_at_detector_node(*n);
    });
}

void GraphFlooder::set_region_frozen(GraphFillRegion &region) {
    bool was_shrinking = region.radius.is_shrinking();
    region.radius = region.radius.then_frozen_at_time(queue.cur_time);

    // No shrinking event can occur while frozen.
    region.shrink_event_tracker.set_no_desired_event();

    // Node events can occur while frozen, from other regions growing into this one.
    // However, those events can only be sooner than the currently scheduled events
    // if the region was previously shrinking (as opposed to growing).
    if (was_shrinking) {
        region.do_op_for_each_node_in_total_area([&](DetectorNode *n) {
            reschedule_events_at_detector_node(*n);
        });
    }
}

void GraphFlooder::set_region_shrinking(GraphFillRegion &region) {
    region.radius = region.radius.then_shrinking_at_time(queue.cur_time);

    // Shrinking events can now occur.
    schedule_tentative_shrink_event(region);

    // No node events can occur while shrinking.
    region.do_op_for_each_node_in_total_area([&](DetectorNode *n) {
        n->node_event_tracker.set_no_desired_event();
    });
}

MwpmEvent GraphFlooder::do_look_at_node_event(DetectorNode &node) {
    auto next = find_next_event_at_node_returning_neighbor_index_and_time(node);
    if (next.second == queue.cur_time) {
        // Need to revisit this node immediately after the mwpm event is handled. There may be an event to handle
        // along another edge, or even along the same edge at a later time (e.g. if this event isn't the *first* event
        // for the neighbor along the current edge when the neighbor is being rescheduled).
        node.node_event_tracker.set_desired_event(
            {
                &node,
                cyclic_time_int{queue.cur_time},
            },
            queue);

        if (node.neighbors[next.first] == nullptr) {
            return do_region_hit_boundary_interaction(node);
        }
        auto &neighbor = *node.neighbors[next.first];
        return do_neighbor_interaction(node, next.first, neighbor);
    } else if (next.first != SIZE_MAX) {
        // Need to revisit this node at a later time.
        node.node_event_tracker.set_desired_event(
            {
                &node,
                cyclic_time_int{next.second},
            },
            queue);
    }

    return MwpmEvent::no_event();
}

MwpmEvent GraphFlooder::process_tentative_event_returning_mwpm_event(FloodCheckEvent tentative_event) {
    switch (tentative_event.tentative_event_type) {
        case LOOK_AT_NODE: {
            return do_look_at_node_event(*tentative_event.data_look_at_node);
        }
        case LOOK_AT_SHRINKING_REGION: {
            return do_region_shrinking(*tentative_event.data_look_at_shrinking_region);
        }
        default:
            throw std::invalid_argument("Unknown tentative event type.");
    }
}

MwpmEvent GraphFlooder::run_until_next_mwpm_notification() {
    while (true) {
        FloodCheckEvent tentative_event = dequeue_valid();
        if (tentative_event.tentative_event_type == NO_FLOOD_CHECK_EVENT) {
            return MwpmEvent::no_event();
        }
        MwpmEvent notification = process_tentative_event_returning_mwpm_event(tentative_event);
        if (notification.event_type != NO_EVENT) {
            return notification;
        }
    }
}

void GraphFlooder::sync_negative_weight_observables_and_detection_events() {
    /// Move set of negative weight detection events into a sorted vector, for faster processing during decoding
    negative_weight_detection_events.clear();
    negative_weight_detection_events.reserve(graph.negative_weight_detection_events_set.size());
    for (auto &det : graph.negative_weight_detection_events_set)
        negative_weight_detection_events.push_back(det);
    /// Move set of negative weight observable indices into a sorted vector, for faster processing during decoding
    negative_weight_observables.clear();
    negative_weight_observables.reserve(graph.negative_weight_observables_set.size());
    for (auto &obs : graph.negative_weight_observables_set)
        negative_weight_observables.push_back(obs);
    negative_weight_obs_mask = 0;
    if (graph.num_observables <= sizeof(pm::obs_int) * 8) {
        // Compute the observable mask for 64 or fewer observables
        for (auto &i : negative_weight_observables)
            negative_weight_obs_mask ^= (pm::obs_int)1 << i;
    }
    negative_weight_sum = graph.negative_weight_sum;
}

GraphFlooder::GraphFlooder() : negative_weight_obs_mask(0), negative_weight_sum(0) {
}
