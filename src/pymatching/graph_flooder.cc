#include "pymatching/graph_flooder.h"
#include "pymatching/graph.h"
#include "pymatching/alternating_tree.h"
#include "pymatching/graph_fill_region.h"
#include "pymatching/varying.h"

pm::GraphFlooder::GraphFlooder(pm::MatchingGraph graph) : graph(std::move(graph)), next_event_vid(0) {
}

pm::GraphFlooder::GraphFlooder(pm::GraphFlooder &&flooder) noexcept
    : graph(std::move(flooder.graph)), queue(std::move(flooder.queue)), next_event_vid(flooder.next_event_vid) {
}

void pm::GraphFlooder::create_region(DetectorNode *node) {
    auto region = new GraphFillRegion();
    auto alt_tree_node = new AltTreeNode(region);
    region->alt_tree_node = alt_tree_node;
    do_region_created_at_empty_detector_node(*region, *node);
}

void pm::GraphFlooder::do_region_created_at_empty_detector_node(GraphFillRegion &region, DetectorNode &detector_node) {
    detector_node.reached_from_source = &detector_node;
    detector_node.distance_from_source = 0;
    detector_node.region_that_arrived = &region;
    region.shell_area.push_back(&detector_node);
    reschedule_events_at_detector_node(detector_node);
}

std::pair<size_t, pm::cumulative_time_int> pm::GraphFlooder::find_next_event_at_node_returning_neighbor_index_and_time(DetectorNode &detector_node) const {
    pm::cumulative_time_int best_time = std::numeric_limits<pm::cumulative_time_int>::max();
    size_t best_neighbor = SIZE_MAX;

    auto rad1 = detector_node.local_radius();

    size_t start = 0;
    if (!detector_node.neighbors.empty() && detector_node.neighbors[0] == nullptr) {
        // If growing towards boundary
        if (rad1.is_growing()) {
            auto weight = detector_node.neighbor_weights[0];
            auto collision_time = (rad1 - weight).time_of_x_intercept_for_growing();
            if (collision_time >= queue.cur_time && collision_time < best_time) {
                best_time = collision_time;
                best_neighbor = 0;
            }
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
        if (!rad1.colliding_with(rad2)) {
            continue;
        }

        auto collision_time = rad1.time_of_x_intercept_when_added_to(rad2 - weight);
        if (collision_time >= queue.cur_time && collision_time < best_time) {
            best_time = collision_time;
            best_neighbor = i;
        }
    }

    return {best_neighbor, best_time};
}

void pm::GraphFlooder::reschedule_events_at_detector_node(DetectorNode &detector_node) {
    auto x = find_next_event_at_node_returning_neighbor_index_and_time(detector_node);
    if (x.first == SIZE_MAX) {
        detector_node.node_event_tracker.set_no_desired_event();
    } else {
        detector_node.node_event_tracker.set_desired_event({
            TentativeEventData_LookAtNode{&detector_node},
            cyclic_time_int{x.second},
        }, queue);
    }
}

void pm::GraphFlooder::reschedule_events_for_region(pm::GraphFillRegion &region) {
    region.shrink_event_tracker.set_no_desired_event();

    if (region.radius.is_shrinking()) {
        schedule_tentative_shrink_event(region);
        region.do_op_for_each_node_in_total_area([](DetectorNode *n) {
            n->invalidate_involved_schedule_items();
        });
    } else {
        region.do_op_for_each_node_in_total_area([this](DetectorNode *n) {
            reschedule_events_at_detector_node(*n);
        });
    }
}

void pm::GraphFlooder::schedule_tentative_shrink_event(pm::GraphFillRegion &region) {
    pm::cumulative_time_int t;
    if (region.shell_area.empty()) {
        t = region.radius.time_of_x_intercept_for_shrinking();
    } else {
        t = region.shell_area.back()->local_radius().time_of_x_intercept_for_shrinking();
    }
    region.shrink_event_tracker.set_desired_event({
        pm::TentativeEventData_LookAtShrinkingRegion{&region},
        cyclic_time_int{t},
    }, queue);
}

void pm::GraphFlooder::do_region_arriving_at_empty_detector_node(
    pm::GraphFillRegion &region, pm::DetectorNode &empty_node, pm::DetectorNode &from_node, size_t neighbor_index) {
    empty_node.observables_crossed_from_source =
        (from_node.observables_crossed_from_source ^ from_node.neighbor_observables[neighbor_index]);
    empty_node.reached_from_source = from_node.reached_from_source;
    empty_node.distance_from_source = from_node.distance_from_source + from_node.neighbor_weights[neighbor_index];
    empty_node.region_that_arrived = &region;
    region.shell_area.push_back(&empty_node);
    reschedule_events_at_detector_node(empty_node);
}

pm::MwpmEvent pm::GraphFlooder::do_region_shrinking(const pm::TentativeEventData_LookAtShrinkingRegion &event) {
    if (event.region->shell_area.empty()) {
        return do_blossom_shattering(*event.region);
    } else if (event.region->shell_area.size() == 1 && event.region->blossom_children.empty()) {
        return do_degenerate_implosion(*event.region);
    } else {
        auto leaving_node = event.region->shell_area.back();
        event.region->shell_area.pop_back();
        leaving_node->region_that_arrived = nullptr;
        leaving_node->reached_from_source = nullptr;
        leaving_node->distance_from_source = -1;
        leaving_node->observables_crossed_from_source = 0;
        reschedule_events_at_detector_node(*leaving_node);
        schedule_tentative_shrink_event(*event.region);
        return MwpmEvent::no_event();
    }
}

pm::MwpmEvent pm::GraphFlooder::do_neighbor_interaction(
    DetectorNode &src,
    size_t src_to_dst_index,
    DetectorNode &dst,
    size_t dst_to_src_index) {
    // First check if one region is moving into an empty location
    if (src.region_that_arrived && !dst.region_that_arrived) {
        do_region_arriving_at_empty_detector_node(
            *src.top_region(),
            dst,
            src,
            src_to_dst_index);
        return MwpmEvent::no_event();
    } else if (dst.region_that_arrived && !src.region_that_arrived) {
        do_region_arriving_at_empty_detector_node(
            *dst.top_region(),
            src,
            dst,
            dst_to_src_index);
        return MwpmEvent::no_event();
    } else {
        // Two regions colliding
        return RegionHitRegionEventData{
            src.top_region(),
            dst.top_region(),
            CompressedEdge(
                src.reached_from_source,
                dst.reached_from_source,
                src.observables_crossed_from_source ^
                    dst.observables_crossed_from_source ^
                    src.neighbor_observables[src_to_dst_index]),
        };
    }
}

pm::MwpmEvent pm::GraphFlooder::do_region_hit_boundary_interaction(DetectorNode &node) {
    return pm::RegionHitBoundaryEventData{
        node.top_region(),
        CompressedEdge(
            node.reached_from_source,
            nullptr,
            node.observables_crossed_from_source ^ node.neighbor_observables[0])};
}

pm::MwpmEvent pm::GraphFlooder::do_degenerate_implosion(const pm::GraphFillRegion &region) {
    return pm::RegionHitRegionEventData{
        region.alt_tree_node->parent.alt_tree_node->outer_region,
        region.alt_tree_node->outer_region,
        pm::CompressedEdge(
            region.alt_tree_node->parent.edge.loc_to,
            region.alt_tree_node->inner_to_outer_edge.loc_to,
            region.alt_tree_node->inner_to_outer_edge.obs_mask ^ region.alt_tree_node->parent.edge.obs_mask)};
}

pm::MwpmEvent pm::GraphFlooder::do_blossom_shattering(pm::GraphFillRegion &region) {
    for (auto &child : region.blossom_children)
        child.region->blossom_parent = nullptr;

    return pm::BlossomShatterEventData{
        &region,
        region.alt_tree_node->parent.edge.loc_from->top_region(),
        region.alt_tree_node->inner_to_outer_edge.loc_from->top_region()};
}

pm::GraphFillRegion *pm::GraphFlooder::create_blossom(std::vector<RegionEdge> &contained_regions) {
    auto blossom_region = new GraphFillRegion();
    blossom_region->radius = pm::Varying32::growing_varying_with_zero_distance_at_time(queue.cur_time);
    blossom_region->blossom_children = std::move(contained_regions);
    for (auto &region_edge : blossom_region->blossom_children) {
        region_edge.region->radius = region_edge.region->radius.then_frozen_at_time(queue.cur_time);
        region_edge.region->blossom_parent = blossom_region;
        region_edge.region->shrink_event_tracker.set_no_desired_event();
    }
    reschedule_events_for_region(*blossom_region);
    return blossom_region;
}

bool pm::GraphFlooder::dequeue_decision(pm::TentativeEvent ev) {
    switch (ev.tentative_event_type) {
        case pm::TentativeType::LOOK_AT_NODE: {
            auto &d = ev.data_look_at_node;
            return d.detector_node->node_event_tracker.dequeue_decision(ev, queue);
        } case pm::TentativeType::LOOK_AT_SHRINKING_REGION: {
            auto &dat = ev.data_look_at_shrinking_region;
            return dat.region->shrink_event_tracker.dequeue_decision(ev, queue);
        } case pm::TentativeType::NO_TENTATIVE_EVENT:
            return true;
        default:
            throw std::invalid_argument("Unrecognized event type.");
    }
}

pm::TentativeEvent pm::GraphFlooder::dequeue_valid() {
    while (true) {
        TentativeEvent ev = queue.dequeue();
        if (dequeue_decision(ev)) {
            return ev;
        }
    }
}

void pm::GraphFlooder::set_region_growing(pm::GraphFillRegion &region) {
    region.radius = region.radius.then_growing_at_time(queue.cur_time);
    reschedule_events_for_region(region);
}

void pm::GraphFlooder::set_region_frozen(pm::GraphFillRegion &region) {
    region.radius = region.radius.then_frozen_at_time(queue.cur_time);
    reschedule_events_for_region(region);
}

void pm::GraphFlooder::set_region_shrinking(pm::GraphFillRegion &region) {
    region.radius = region.radius.then_shrinking_at_time(queue.cur_time);
    reschedule_events_for_region(region);
}


pm::MwpmEvent pm::GraphFlooder::do_valid_tentative_event_returning_mwpm_event(TentativeEvent tentative_event) {
    switch (tentative_event.tentative_event_type) {
        case LOOK_AT_NODE: {
            auto &node = *tentative_event.data_look_at_node.detector_node;

            auto ev = find_next_event_at_node_returning_neighbor_index_and_time(node);
            if (ev.second != queue.cur_time) {
                // The event was invalidated by a change between enqueuing and now.
                return pm::MwpmEvent::no_event();
            }

            auto *neighbor = node.neighbors[ev.first];
            if (neighbor == nullptr) {
                return do_region_hit_boundary_interaction(node);
            } else {
                size_t back_index = SIZE_MAX >> 2;  // Well into segfault land, please.
                for (size_t k = 0; k < neighbor->neighbors.size(); k++) {
                    if (neighbor->neighbors[k] == &node) {
                        back_index = k;
                        break;
                    }
                }
                return do_neighbor_interaction(node, ev.first, *neighbor, back_index);
            }
        } case LOOK_AT_SHRINKING_REGION: {
            return do_region_shrinking(tentative_event.data_look_at_shrinking_region);
        } default:
            throw std::invalid_argument("Unknown tentative event type.");
    }
}

pm::MwpmEvent pm::GraphFlooder::next_event() {
    while (true) {
        TentativeEvent tentative_event = dequeue_valid();
        if (tentative_event.tentative_event_type == NO_TENTATIVE_EVENT) {
            return MwpmEvent::no_event();
        }
        MwpmEvent processed = do_valid_tentative_event_returning_mwpm_event(tentative_event);
        if (tentative_event.tentative_event_type == LOOK_AT_NODE) {
            reschedule_events_at_detector_node(*tentative_event.data_look_at_node.detector_node);
        }
        if (processed.event_type != NO_EVENT) {
            return processed;
        }
    }
}
