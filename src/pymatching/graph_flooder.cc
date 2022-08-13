#include "pymatching/graph_flooder.h"

#include "pymatching/graph.h"
#include "pymatching/graph_fill_region.h"
#include "pymatching/varying.h"

pm::GraphFlooder::GraphFlooder(pm::MatchingGraph graph, size_t num_buckets) : graph(std::move(graph)), queue(num_buckets), next_event_vid(0) {
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

void pm::GraphFlooder::reschedule_events_at_detector_node(DetectorNode &detector_node) {
    detector_node.invalidate_involved_schedule_items();
    auto rad1 = detector_node.local_radius();
    for (size_t i = 0; i < detector_node.neighbors.size(); i++) {
        auto weight = detector_node.neighbor_weights[i];
        if (!detector_node.neighbors[i]) {
            // If growing towards boundary
            if (rad1.is_growing()) {
                auto collision_time = (rad1 - weight).time_of_x_intercept_for_growing();
                schedule_tentative_neighbor_interaction_event(
                    &detector_node, i, nullptr, -1, collision_time);
            }
            continue;
        }
        auto neighbor = detector_node.neighbors[i];
        if (detector_node.has_same_owner_as(*neighbor))
            continue;
        auto rad2 = neighbor->local_radius();
        if (!rad1.colliding_with(rad2))
            continue;
        auto collision_time = rad1.time_of_x_intercept_when_added_to(rad2 - weight);
        if (collision_time < (int64_t)queue.cur_time) {
            continue;
        }

        // Find index of detector_node from neighbor
        size_t j = 0;
        for (auto &n : neighbor->neighbors) {
            if (n == &detector_node)
                break;
            j++;
        }
        schedule_tentative_neighbor_interaction_event(&detector_node, i, neighbor, j, collision_time);
    }
}

void pm::GraphFlooder::reschedule_events_for_region(pm::GraphFillRegion &region) {
    // Invalidate existing events by setting the vid to an index that doesn't correspond to an event
    // affecting the object. Note that decrementing would not be safe due to ABA writes.
    region.shrink_event_vid++;

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

void pm::GraphFlooder::schedule_tentative_neighbor_interaction_event(
    pm::DetectorNode *detector_node_1,
    size_t schedule_list_index_1,
    pm::DetectorNode *detector_node_2,
    size_t schedule_list_index_2,
    pm::time_int event_time) {
    auto vid = next_event_vid++;
    detector_node_1->edge_event_vids[schedule_list_index_1] = vid;
    if (detector_node_2 != nullptr) {
        detector_node_2->edge_event_vids[schedule_list_index_2] = vid;
    }
    queue.enqueue({
        pm::TentativeNeighborInteractionEventData{
            detector_node_1, schedule_list_index_1, detector_node_2, schedule_list_index_2
        },
        event_time,
        vid,
    });
}

void pm::GraphFlooder::schedule_tentative_shrink_event(pm::GraphFillRegion &region) {
    pm::time_int t;
    if (region.shell_area.empty()) {
        t = region.radius.time_of_x_intercept_for_shrinking();
    } else {
        t = region.shell_area.back()->local_radius().time_of_x_intercept_for_shrinking();
    }
    auto vid = next_event_vid++;
    region.shrink_event_vid = vid;
    queue.enqueue({
        pm::TentativeRegionShrinkEventData{
            &region,
        },
        t,
        vid,
    });
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

pm::MwpmEvent pm::GraphFlooder::do_region_shrinking(const pm::TentativeRegionShrinkEventData &event) {
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

pm::MwpmEvent pm::GraphFlooder::do_neighbor_interaction(const pm::TentativeNeighborInteractionEventData &event) {
    // First check if one region is moving into an empty location
    if (event.detector_node_1->region_that_arrived && !event.detector_node_2->region_that_arrived) {
        do_region_arriving_at_empty_detector_node(
            *event.detector_node_1->top_region(),
            *event.detector_node_2,
            *event.detector_node_1,
            event.node_1_neighbor_index);
        return MwpmEvent::no_event();
    } else if (event.detector_node_2->region_that_arrived && !event.detector_node_1->region_that_arrived) {
        do_region_arriving_at_empty_detector_node(
            *event.detector_node_2->top_region(),
            *event.detector_node_1,
            *event.detector_node_2,
            event.node_2_neighbor_index);
        return MwpmEvent::no_event();
    } else {
        // Two regions colliding
        return RegionHitRegionEventData{
            event.detector_node_1->top_region(),
            event.detector_node_2->top_region(),
            CompressedEdge(
                event.detector_node_1->reached_from_source,
                event.detector_node_2->reached_from_source,
                event.detector_node_1->observables_crossed_from_source ^
                    event.detector_node_2->observables_crossed_from_source ^
                    event.detector_node_1->neighbor_observables[event.node_1_neighbor_index]),
        };
    }
}

pm::MwpmEvent pm::GraphFlooder::do_region_hit_boundary_interaction(
    const pm::TentativeNeighborInteractionEventData &event) {
    return pm::RegionHitBoundaryEventData{
        event.detector_node_1->top_region(),
        CompressedEdge(
            event.detector_node_1->reached_from_source,
            nullptr,
            event.detector_node_1->observables_crossed_from_source ^
                event.detector_node_1->neighbor_observables[event.node_1_neighbor_index])};
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
        region_edge.region->shrink_event_vid++; // Invalidate events affecting the region.
    }
    reschedule_events_for_region(*blossom_region);
    return blossom_region;
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
        case INTERACTION:
            if (tentative_event.neighbor_interaction_event_data.detector_node_2) {
                return do_neighbor_interaction(tentative_event.neighbor_interaction_event_data);
            } else {
                return do_region_hit_boundary_interaction(tentative_event.neighbor_interaction_event_data);
            }
        case SHRINKING:
            return do_region_shrinking(tentative_event.region_shrink_event_data);
        default:
            throw std::invalid_argument("Unknown tentative event type.");
    }
}

pm::MwpmEvent pm::GraphFlooder::next_event() {
    while (true) {
        TentativeEvent tentative_event; // NOLINT(cppcoreguidelines-pro-type-member-init)
        if (!queue.try_pop(&tentative_event)) {
            return MwpmEvent::no_event();
        }
        MwpmEvent processed = do_valid_tentative_event_returning_mwpm_event(tentative_event);
        if (processed.event_type != NO_EVENT) {
            return processed;
        }
    }
}
