#include "events.h"


pm::TentativeEvent::TentativeEvent(pm::DetectorNode *detector_node_1, size_t node_1_neighbor_index,
                               pm::DetectorNode *detector_node_2, size_t node_2_neighbor_index, time_int time)
                               : neighbor_interaction_event_data(detector_node_1, node_1_neighbor_index,
                                                                 detector_node_2, node_2_neighbor_index),
                               time(time), tentative_event_type(INTERACTION), is_invalidated(false){}


pm::TentativeEvent::TentativeEvent(pm::GraphFillRegion* region, time_int time)
        : region_shrink_event_data(region),
          time(time), tentative_event_type(SHRINKING), is_invalidated(false){}


void pm::TentativeEvent::invalidate() {
    is_invalidated = true;
    // Is resetting schedule pointers below needed?
    if (tentative_event_type == INTERACTION) {
        neighbor_interaction_event_data.detector_node_1
            ->neighbor_schedules[neighbor_interaction_event_data.node_1_neighbor_index] = nullptr;
        if (neighbor_interaction_event_data.detector_node_2){
            neighbor_interaction_event_data.detector_node_2
                ->neighbor_schedules[neighbor_interaction_event_data.node_2_neighbor_index] = nullptr;
        }
    } else if (tentative_event_type == SHRINKING) {
        region_shrink_event_data.region->shrink_event = nullptr;
    }
}

pm::TentativeNeighborInteractionEventData::TentativeNeighborInteractionEventData(pm::DetectorNode *detector_node_1,
                                                                                 size_t node_1_neighbor_index,
                                                                                 pm::DetectorNode *detector_node_2,
                                                                                 size_t node_2_neighbor_index)
                                           : detector_node_1(detector_node_1),
                                           node_1_neighbor_index(node_1_neighbor_index),
                                           detector_node_2(detector_node_2),
                                           node_2_neighbor_index(node_2_neighbor_index){}

pm::TentativeRegionShrinkEventData::TentativeRegionShrinkEventData(pm::GraphFillRegion* region)
    : region(region) {}
