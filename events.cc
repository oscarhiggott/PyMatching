#include "events.h"


pm::TentativeEvent::TentativeEvent(pm::DetectorNode *detector_node_1, size_t node_1_neighbor_index,
                               pm::DetectorNode *detector_node_2, size_t node_2_neighbor_index, time_int time, uint32_t validation_index)
                               : neighbor_interaction_event_data(detector_node_1, node_1_neighbor_index,
                                                                 detector_node_2, node_2_neighbor_index),
                               time(time), tentative_event_type(INTERACTION), validation_index(validation_index){}


pm::TentativeEvent::TentativeEvent(pm::GraphFillRegion* region, time_int time, uint32_t validation_index)
        : region_shrink_event_data(region),
          time(time), tentative_event_type(SHRINKING), validation_index(validation_index){}


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
        validation_index != rhs.validation_index)
        return false;
    switch (tentative_event_type) {
        case SHRINKING:
            return region_shrink_event_data == rhs.region_shrink_event_data;
        case INTERACTION:
            return neighbor_interaction_event_data == rhs.neighbor_interaction_event_data;
        default:
            throw std::invalid_argument("Unrecognized event type");
    }
}

bool pm::TentativeEvent::operator!=(const TentativeEvent &rhs) const {
    return !(rhs == *this);
}

pm::RegionHitRegionEventData::RegionHitRegionEventData(
        pm::GraphFillRegion *region1, pm::GraphFillRegion *region2, pm::CompressedEdge edge
        )
            : region1(region1), region2(region2), edge(edge) {}

bool pm::RegionHitRegionEventData::operator==(const pm::RegionHitRegionEventData &rhs) const {
    return region1 == rhs.region1 &&
           region2 == rhs.region2 &&
           edge == rhs.edge;
}

bool pm::RegionHitRegionEventData::operator!=(const pm::RegionHitRegionEventData &rhs) const {
    return !(rhs == *this);
}

pm::RegionHitBoundaryEventData::RegionHitBoundaryEventData(
        pm::GraphFillRegion *region, CompressedEdge edge
        )
        : region(region), edge(edge) {}

bool pm::RegionHitBoundaryEventData::operator==(const pm::RegionHitBoundaryEventData &rhs) const {
    return region == rhs.region &&
           edge == rhs.edge;
}

bool pm::RegionHitBoundaryEventData::operator!=(const pm::RegionHitBoundaryEventData &rhs) const {
    return !(rhs == *this);
}


pm::BlossomShatterEventData::BlossomShatterEventData(
        pm::GraphFillRegion *blossom_region, pm::GraphFillRegion *in_parent_region,
        pm::GraphFillRegion *in_child_region)
        : blossom_region(blossom_region), in_parent_region(in_parent_region), in_child_region(in_child_region) {}

bool pm::BlossomShatterEventData::operator==(const pm::BlossomShatterEventData &rhs) const {
    return blossom_region == rhs.blossom_region &&
           in_parent_region == rhs.in_parent_region &&
           in_child_region == rhs.in_child_region;
}

bool pm::BlossomShatterEventData::operator!=(const pm::BlossomShatterEventData &rhs) const {
    return !(rhs == *this);
}

pm::MwpmEvent::MwpmEvent(pm::RegionHitRegionEventData region_hit_region_event_data)
    : region_hit_region_event_data(region_hit_region_event_data),
      event_type(pm::REGION_HIT_REGION) {
}
pm::MwpmEvent::MwpmEvent(pm::RegionHitBoundaryEventData region_hit_region_event_data)
    : region_hit_boundary_event_data(region_hit_region_event_data),
      event_type(pm::REGION_HIT_BOUNDARY) {
}
pm::MwpmEvent::MwpmEvent(pm::BlossomShatterEventData region_hit_region_event_data)
    : blossom_shatter_event_data(region_hit_region_event_data),
      event_type(pm::BLOSSOM_SHATTER) {
}
pm::MwpmEvent::MwpmEvent() : event_type(NO_EVENT) {
}

bool pm::MwpmEvent::operator==(const pm::MwpmEvent &rhs) const {
    if (event_type != rhs.event_type)
        return false;
    switch (event_type) {
        case NO_EVENT:
            return true;
        case REGION_HIT_REGION:
            return region_hit_region_event_data == rhs.region_hit_region_event_data;
        case REGION_HIT_BOUNDARY:
            return region_hit_boundary_event_data == rhs.region_hit_boundary_event_data;
        case BLOSSOM_SHATTER:
            return blossom_shatter_event_data == rhs.blossom_shatter_event_data;
        default:
            throw std::invalid_argument("Unrecognized event type");
    }
}

bool pm::MwpmEvent::operator!=(const pm::MwpmEvent &rhs) const {
    return !(rhs == *this);
}
