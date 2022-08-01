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

bool pm::TentativeNeighborInteractionEventData::operator==(const pm::TentativeNeighborInteractionEventData &rhs) const {
    return detector_node_1 == rhs.detector_node_1 &&
           node_1_neighbor_index == rhs.node_1_neighbor_index &&
           detector_node_2 == rhs.detector_node_2 &&
           node_2_neighbor_index == rhs.node_2_neighbor_index;
}

bool pm::TentativeNeighborInteractionEventData::operator!=(const pm::TentativeNeighborInteractionEventData &rhs) const {
    return !(rhs == *this);
}

pm::TentativeRegionShrinkEventData::TentativeRegionShrinkEventData(pm::GraphFillRegion* region)
    : region(region) {}

bool pm::TentativeRegionShrinkEventData::operator==(const pm::TentativeRegionShrinkEventData &rhs) const {
    return region == rhs.region;
}

bool pm::TentativeRegionShrinkEventData::operator!=(const pm::TentativeRegionShrinkEventData &rhs) const {
    return !(rhs == *this);
}

bool pm::TentativeEvent::operator==(const TentativeEvent &rhs) const {
    if (time != rhs.time || tentative_event_type != rhs.tentative_event_type ||
        is_invalidated != rhs.is_invalidated)
        return false;
    switch (tentative_event_type) {
        case SHRINKING:
            return region_shrink_event_data == rhs.region_shrink_event_data;
        case INTERACTION:
            return neighbor_interaction_event_data == rhs.neighbor_interaction_event_data;
    }
}

bool pm::TentativeEvent::operator!=(const TentativeEvent &rhs) const {
    return !(rhs == *this);
}