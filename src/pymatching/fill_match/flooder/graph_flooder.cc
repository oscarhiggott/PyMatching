#include "pymatching/fill_match/flooder/graph_flooder.h"

#include "pymatching/fill_match/flooder/graph.h"
#include "pymatching/fill_match/flooder/graph_fill_region.h"
#include "pymatching/fill_match/flooder/varying.h"
#include "pymatching/fill_match/matcher/alternating_tree.h"

using namespace pm;

GraphFlooder::GraphFlooder(MatchingGraph graph) : graph(std::move(graph)), next_event_vid(0) {
}

GraphFlooder::GraphFlooder(GraphFlooder &&flooder) noexcept
    : graph(std::move(flooder.graph)), queue(std::move(flooder.queue)), next_event_vid(flooder.next_event_vid) {
}

void GraphFlooder::create_region(DetectorNode *node) {
    auto region = new GraphFillRegion();
    auto alt_tree_node = new AltTreeNode(region);
    region->alt_tree_node = alt_tree_node;
    do_region_created_at_empty_detector_node(*region, *node);
}

void GraphFlooder::do_region_created_at_empty_detector_node(GraphFillRegion &region, DetectorNode &detector_node) {
    detector_node.reached_from_source = &detector_node;
    detector_node.distance_from_source = 0;
    detector_node.region_that_arrived = &region;
    region.shell_area.push_back(&detector_node);
    reschedule_events_at_detector_node(detector_node);
}

std::pair<size_t, cumulative_time_int> GraphFlooder::find_next_event_at_node_returning_neighbor_index_and_time(
    DetectorNode &detector_node) const {
    cumulative_time_int best_time = std::numeric_limits<cumulative_time_int>::max();
    size_t best_neighbor = SIZE_MAX;

    //    auto rad1 = detector_node.total_radius();
    //    if (rad1.is_shrinking()) {
    //        // No node collision events can occur while shrinking.
    //        return {best_neighbor, best_time};
    //    }
    //
    //    rad1 -= detector_node.distance_from_source;
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

void GraphFlooder::reschedule_events_at_detector_node(DetectorNode &detector_node) {
    auto x = find_next_event_at_node_returning_neighbor_index_and_time(detector_node);
    if (x.first == SIZE_MAX) {
        detector_node.node_event_tracker.set_no_desired_event();
    } else {
        detector_node.node_event_tracker.set_desired_event(
            {
                TentativeEventData_LookAtNode{&detector_node},
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
            TentativeEventData_LookAtShrinkingRegion{&region},
            cyclic_time_int{t},
        },
        queue);
}

void GraphFlooder::do_region_arriving_at_empty_detector_node(
    GraphFillRegion &region, DetectorNode &empty_node, DetectorNode &from_node, size_t neighbor_index) {
    empty_node.observables_crossed_from_source =
        (from_node.observables_crossed_from_source ^ from_node.neighbor_observables[neighbor_index]);
    empty_node.reached_from_source = from_node.reached_from_source;
    empty_node.distance_from_source = from_node.distance_from_source + from_node.neighbor_weights[neighbor_index];
    empty_node.region_that_arrived = &region;
    region.shell_area.push_back(&empty_node);
    reschedule_events_at_detector_node(empty_node);
}

MwpmEvent GraphFlooder::do_region_shrinking(const TentativeEventData_LookAtShrinkingRegion &event) {
    if (event.region->shell_area.empty()) {
        return do_blossom_shattering(*event.region);
    } else if (event.region->shell_area.size() == 1 && event.region->blossom_children.empty()) {
        return do_degenerate_implosion(*event.region);
    } else {
        auto leaving_node = event.region->shell_area.back();
        event.region->shell_area.pop_back();
        leaving_node->region_that_arrived = nullptr;
        leaving_node->reached_from_source = nullptr;
        leaving_node->distance_from_source = 0;
        leaving_node->observables_crossed_from_source = 0;
        reschedule_events_at_detector_node(*leaving_node);
        schedule_tentative_shrink_event(*event.region);
        return MwpmEvent::no_event();
    }
}

MwpmEvent GraphFlooder::do_neighbor_interaction(
    DetectorNode &src, size_t src_to_dst_index, DetectorNode &dst, size_t dst_to_src_index) {
    // First check if one region is moving into an empty location
    if (src.region_that_arrived && !dst.region_that_arrived) {
        do_region_arriving_at_empty_detector_node(*src.top_region(), dst, src, src_to_dst_index);
        return MwpmEvent::no_event();
    } else if (dst.region_that_arrived && !src.region_that_arrived) {
        do_region_arriving_at_empty_detector_node(*dst.top_region(), src, dst, dst_to_src_index);
        return MwpmEvent::no_event();
    } else {
        // Two regions colliding
        return RegionHitRegionEventData{
            src.top_region(),
            dst.top_region(),
            CompressedEdge(
                src.reached_from_source,
                dst.reached_from_source,
                src.observables_crossed_from_source ^ dst.observables_crossed_from_source ^
                    src.neighbor_observables[src_to_dst_index]),
        };
    }
}

MwpmEvent GraphFlooder::do_region_hit_boundary_interaction(DetectorNode &node) {
    return RegionHitBoundaryEventData{
        node.top_region(),
        CompressedEdge(
            node.reached_from_source, nullptr, node.observables_crossed_from_source ^ node.neighbor_observables[0])};
}

MwpmEvent GraphFlooder::do_degenerate_implosion(const GraphFillRegion &region) {
    return RegionHitRegionEventData{
        region.alt_tree_node->parent.alt_tree_node->outer_region,
        region.alt_tree_node->outer_region,
        CompressedEdge(
            region.alt_tree_node->parent.edge.loc_to,
            region.alt_tree_node->inner_to_outer_edge.loc_to,
            region.alt_tree_node->inner_to_outer_edge.obs_mask ^ region.alt_tree_node->parent.edge.obs_mask)};
}

MwpmEvent GraphFlooder::do_blossom_shattering(GraphFillRegion &region) {
    for (auto &child : region.blossom_children)
        child.region->blossom_parent = nullptr;

    return BlossomShatterEventData{
        &region,
        region.alt_tree_node->parent.edge.loc_from->top_region(),
        region.alt_tree_node->inner_to_outer_edge.loc_from->top_region()};
}

GraphFillRegion *GraphFlooder::create_blossom(std::vector<RegionEdge> &contained_regions) {
    auto blossom_region = new GraphFillRegion();
    blossom_region->radius = Varying32::growing_varying_with_zero_distance_at_time(queue.cur_time);
    blossom_region->blossom_children = std::move(contained_regions);
    for (auto &region_edge : blossom_region->blossom_children) {
        region_edge.region->radius = region_edge.region->radius.then_frozen_at_time(queue.cur_time);
        region_edge.region->blossom_parent = blossom_region;
        region_edge.region->shrink_event_tracker.set_no_desired_event();
    }

    blossom_region->do_op_for_each_node_in_total_area([this](DetectorNode *n) {
        reschedule_events_at_detector_node(*n);
    });

    return blossom_region;
}

bool GraphFlooder::dequeue_decision(FloodCheckEvent ev) {
    switch (ev.tentative_event_type) {
        case TentativeType::LOOK_AT_NODE: {
            auto &d = ev.data_look_at_node;
            return d.detector_node->node_event_tracker.dequeue_decision(ev, queue);
        }
        case TentativeType::LOOK_AT_SHRINKING_REGION: {
            auto &dat = ev.data_look_at_shrinking_region;
            return dat.region->shrink_event_tracker.dequeue_decision(ev, queue);
        }
        case TentativeType::NO_TENTATIVE_EVENT:
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
        // Need to revisit this node immediately after the mwpm event is handled.
        node.node_event_tracker.set_desired_event(
            {
                TentativeEventData_LookAtNode{&node},
                cyclic_time_int{queue.cur_time},
            },
            queue);

        if (node.neighbors[next.first] == nullptr) {
            return do_region_hit_boundary_interaction(node);
        }
        auto &neighbor = *node.neighbors[next.first];
        return do_neighbor_interaction(node, next.first, neighbor, neighbor.index_of_neighbor(&node));
    } else if (next.first != SIZE_MAX) {
        // Need to revisit this node at a later time.
        node.node_event_tracker.set_desired_event(
            {
                TentativeEventData_LookAtNode{&node},
                cyclic_time_int{next.second},
            },
            queue);
    }

    return MwpmEvent::no_event();
}

MwpmEvent GraphFlooder::process_tentative_event_returning_mwpm_event(FloodCheckEvent tentative_event) {
    switch (tentative_event.tentative_event_type) {
        case LOOK_AT_NODE: {
            return do_look_at_node_event(*tentative_event.data_look_at_node.detector_node);
        }
        case LOOK_AT_SHRINKING_REGION: {
            return do_region_shrinking(tentative_event.data_look_at_shrinking_region);
        }
        default:
            throw std::invalid_argument("Unknown tentative event type.");
    }
}

MwpmEvent GraphFlooder::run_until_next_mwpm_notification() {
    while (true) {
        FloodCheckEvent tentative_event = dequeue_valid();
        if (tentative_event.tentative_event_type == NO_TENTATIVE_EVENT) {
            return MwpmEvent::no_event();
        }
        MwpmEvent notification = process_tentative_event_returning_mwpm_event(tentative_event);
        if (notification.event_type != NO_EVENT) {
            return notification;
        }
    }
}
