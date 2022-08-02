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


pm::BlossomImplodeEventData::BlossomImplodeEventData(
        pm::GraphFillRegion *blossom_region, pm::GraphFillRegion *in_parent_region,
        pm::GraphFillRegion *in_child_region)
        : blossom_region(blossom_region), in_parent_region(in_parent_region), in_child_region(in_child_region) {}

bool pm::BlossomImplodeEventData::operator==(const pm::BlossomImplodeEventData &rhs) const {
    return blossom_region == rhs.blossom_region &&
           in_parent_region == rhs.in_parent_region &&
           in_child_region == rhs.in_child_region;
}

bool pm::BlossomImplodeEventData::operator!=(const pm::BlossomImplodeEventData &rhs) const {
    return !(rhs == *this);
}

pm::MwpmEvent::MwpmEvent(pm::GraphFillRegion *region1, pm::GraphFillRegion *region2, pm::CompressedEdge edge)
    : region_hit_region_event_data(region1, region2, edge), event_type(REGION_HIT_REGION) {}

pm::MwpmEvent::MwpmEvent(pm::GraphFillRegion *region, pm::CompressedEdge edge)
    :  region_hit_boundary_event_data(region, edge), event_type(REGION_HIT_BOUNDARY) {}

pm::MwpmEvent::MwpmEvent(pm::GraphFillRegion *blossom_region, pm::GraphFillRegion *in_parent_region,
                         pm::GraphFillRegion *in_child_region)
     : blossom_implode_event_data(blossom_region, in_parent_region, in_child_region), event_type(BLOSSOM_IMPLODE){}

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
        case BLOSSOM_IMPLODE:
            return blossom_implode_event_data == rhs.blossom_implode_event_data;
    }

    return region_hit_region_event_data == rhs.region_hit_region_event_data &&
           region_hit_boundary_event_data == rhs.region_hit_boundary_event_data &&
           blossom_implode_event_data == rhs.blossom_implode_event_data &&
           event_type == rhs.event_type;
}

bool pm::MwpmEvent::operator!=(const pm::MwpmEvent &rhs) const {
    return !(rhs == *this);
}
