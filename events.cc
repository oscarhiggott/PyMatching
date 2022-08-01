#include "events.h"

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
        tentative_region_shrink_event_data.region->shrink_event = nullptr;
    }
}

