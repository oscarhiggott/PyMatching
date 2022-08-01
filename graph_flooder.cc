#include "graph.h"
#include "graph_flooder.h"
#include "graph_fill_region.h"


pm::GraphFlooder::GraphFlooder(pm::Graph &graph) : graph(std::move(graph)), time(0) {};

void pm::GraphFlooder::create_region(DetectorNode *node) {
    auto region = new GraphFillRegion();
    auto alt_tree_node =  new AltTreeNode(region);
    region->alt_tree_node = alt_tree_node;
    do_region_created_at_empty_detector_node(*region, *node);
}

void pm::GraphFlooder::do_region_created_at_empty_detector_node(GraphFillRegion &region,
                                                                DetectorNode &detector_node) {
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
                schedule_tentative_neighbor_interaction_event(
                        &detector_node,
                        i,
                        nullptr,
                        -1,
                        (rad1 - weight).time_of_x_intercept_for_growing()
                        );
            }
            continue;
        }
        auto neighbor = detector_node.neighbors[i];
        if (detector_node.has_same_owner_as(*neighbor))
            continue;
        auto rad2 = neighbor->local_radius();
        if (!rad1.colliding_with(rad2))
            continue;
        auto coll_time = rad1.time_of_x_intercept_when_added_to((rad2 - weight));

        // Find index of detector_node from neighbor
        size_t j = 0;
        for (auto n : neighbor->neighbors){
            if (n == &detector_node)
                break;
            j++;
        }
        schedule_tentative_neighbor_interaction_event(
                &detector_node,
                i,
                neighbor,
                j,
                coll_time
                );
    }
}

void pm::GraphFlooder::schedule_tentative_neighbor_interaction_event(pm::DetectorNode *detector_node_1,
                                                                     size_t schedule_list_index_1,
                                                                     pm::DetectorNode *detector_node_2,
                                                                     size_t schedule_list_index_2,
                                                                     pm::time_int event_time) {
    auto e = new pm::TentativeEvent(detector_node_1, schedule_list_index_1, detector_node_2,
                                    schedule_list_index_2, event_time);
    detector_node_1->neighbor_schedules[schedule_list_index_1] = e;
    if (detector_node_2)
        detector_node_2->neighbor_schedules[schedule_list_index_2] = e;
    queue.push(e);
}
