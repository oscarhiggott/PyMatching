#include "pymatching/events.h"

pm::TentativeEvent::TentativeEvent(
    pm::TentativeNeighborInteractionEventData data,
    time_int time,
    uint64_t validation_index)
    : neighbor_interaction_event_data(data),
      time(time),
      tentative_event_type(INTERACTION),
      vid(validation_index) {
}

pm::TentativeEvent::TentativeEvent(pm::TentativeRegionShrinkEventData data, time_int time, uint64_t validation_index)
    : region_shrink_event_data(data), time(time), tentative_event_type(SHRINKING), vid(validation_index) {
}
pm::TentativeEvent::TentativeEvent(time_int time, uint64_t validation_index) : time(time), tentative_event_type(NO_TENTATIVE_EVENT), vid(validation_index) {
}


bool pm::TentativeEvent::is_still_valid() const {
    switch (tentative_event_type) {
        case INTERACTION: {
            auto &d = neighbor_interaction_event_data;
            if (d.detector_node_1->edge_event_vids[d.node_1_neighbor_index] != vid) {
                return false;
            }
            if (d.detector_node_2 != nullptr && d.detector_node_2->edge_event_vids[d.node_2_neighbor_index] != vid) {
                return false;
            }
            return true;
        } case SHRINKING:
            return region_shrink_event_data.region->shrink_event_vid == vid;
        case NO_TENTATIVE_EVENT:
            return vid == 0;
        default:
            throw std::invalid_argument("Unrecognized event type.");
    }
}

pm::TentativeNeighborInteractionEventData::TentativeNeighborInteractionEventData(
    pm::DetectorNode *detector_node_1,
    size_t node_1_neighbor_index,
    pm::DetectorNode *detector_node_2,
    size_t node_2_neighbor_index)
    : detector_node_1(detector_node_1),
      node_1_neighbor_index(node_1_neighbor_index),
      detector_node_2(detector_node_2),
      node_2_neighbor_index(node_2_neighbor_index) {
}

bool pm::TentativeNeighborInteractionEventData::operator==(const pm::TentativeNeighborInteractionEventData &rhs) const {
    return detector_node_1 == rhs.detector_node_1 && node_1_neighbor_index == rhs.node_1_neighbor_index &&
           detector_node_2 == rhs.detector_node_2 && node_2_neighbor_index == rhs.node_2_neighbor_index;
}

bool pm::TentativeNeighborInteractionEventData::operator!=(const pm::TentativeNeighborInteractionEventData &rhs) const {
    return !(rhs == *this);
}

pm::TentativeRegionShrinkEventData::TentativeRegionShrinkEventData(pm::GraphFillRegion *region) : region(region) {
}

bool pm::TentativeRegionShrinkEventData::operator==(const pm::TentativeRegionShrinkEventData &rhs) const {
    return region == rhs.region;
}

bool pm::TentativeRegionShrinkEventData::operator!=(const pm::TentativeRegionShrinkEventData &rhs) const {
    return !(rhs == *this);
}

bool pm::TentativeEvent::operator==(const TentativeEvent &rhs) const {
    if (time != rhs.time || tentative_event_type != rhs.tentative_event_type || vid != rhs.vid)
        return false;
    switch (tentative_event_type) {
        case SHRINKING:
            return region_shrink_event_data == rhs.region_shrink_event_data;
        case INTERACTION:
            return neighbor_interaction_event_data == rhs.neighbor_interaction_event_data;
        case NO_TENTATIVE_EVENT:
            return true;
        default:
            throw std::invalid_argument("Unrecognized event type");
    }
}

std::ostream &pm::operator<<(std::ostream &out, const TentativeEvent &ev) {
    out << "TentativeEvent{.time=";
    out << ev.time;
    out << ", .vid=";
    out << ev.vid;
    out << ", .type=";
    switch (ev.tentative_event_type) {
        case SHRINKING:
            out << "SHRINKING";
            break;
        case INTERACTION:
            out << "INTERACTION";
            break;
        case NO_TENTATIVE_EVENT:
            out << "NO_TENTATIVE_EVENT";
            break;
        default:
            throw std::invalid_argument("Unrecognized event type");
    }
    out << "}";
    return out;
}

std::ostream &pm::operator<<(std::ostream &out, const pm::RegionHitRegionEventData &dat) {
    out << "{.region1=" << dat.region1;
    out << ", .region2=" << dat.region2;
    out << ", .edge=" << dat.edge;
    out << "}";
    return out;
}

std::ostream &pm::operator<<(std::ostream &out, const pm::BlossomShatterEventData &dat) {
    out << "{.blossom=" << dat.blossom_region;
    out << ", .in_parent_region=" << dat.in_parent_region;
    out << ", .in_child_region=" << dat.in_child_region;
    out << "}";
    return out;
}

std::ostream &pm::operator<<(std::ostream &out, const MwpmEvent &ev) {
    out << "MwpmEvent{.type=";
    switch (ev.event_type) {
        case NO_EVENT:
            out << "NO_EVENT";
            break;
        case REGION_HIT_REGION:
            out << "REGION_HIT_REGION, .dat=" << ev.region_hit_region_event_data;
            break;
        case BLOSSOM_SHATTER:
            out << "BLOSSOM_SHATTER, .dat=" << ev.blossom_shatter_event_data;
            break;
        case REGION_HIT_BOUNDARY:
            out << "REGION_HIT_BOUNDARY";
            break;
        default:
            throw std::invalid_argument("Unrecognized event type");
    }
    out << "}";
    return out;
}

std::string pm::RegionHitRegionEventData::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

std::string pm::BlossomShatterEventData::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

std::string pm::TentativeEvent::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

std::string pm::MwpmEvent::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

bool pm::TentativeEvent::operator!=(const TentativeEvent &rhs) const {
    return !(rhs == *this);
}

pm::RegionHitRegionEventData::RegionHitRegionEventData(
    pm::GraphFillRegion *region1, pm::GraphFillRegion *region2, pm::CompressedEdge edge)
    : region1(region1), region2(region2), edge(edge) {
}

bool pm::RegionHitRegionEventData::operator==(const pm::RegionHitRegionEventData &rhs) const {
    if (region1 == rhs.region1 && region2 == rhs.region2 && edge == rhs.edge) {
        return true;
    }
    // Reversed also okay.
    return region1 == rhs.region2 && region1 == rhs.region2 && edge == rhs.edge.reversed();
}

bool pm::RegionHitRegionEventData::operator!=(const pm::RegionHitRegionEventData &rhs) const {
    return !(rhs == *this);
}

pm::RegionHitBoundaryEventData::RegionHitBoundaryEventData(pm::GraphFillRegion *region, CompressedEdge edge)
    : region(region), edge(edge) {
}

bool pm::RegionHitBoundaryEventData::operator==(const pm::RegionHitBoundaryEventData &rhs) const {
    return region == rhs.region && edge == rhs.edge;
}

bool pm::RegionHitBoundaryEventData::operator!=(const pm::RegionHitBoundaryEventData &rhs) const {
    return !(rhs == *this);
}

pm::BlossomShatterEventData::BlossomShatterEventData(
    pm::GraphFillRegion *blossom_region, pm::GraphFillRegion *in_parent_region, pm::GraphFillRegion *in_child_region)
    : blossom_region(blossom_region), in_parent_region(in_parent_region), in_child_region(in_child_region) {
}

bool pm::BlossomShatterEventData::operator==(const pm::BlossomShatterEventData &rhs) const {
    return blossom_region == rhs.blossom_region && in_parent_region == rhs.in_parent_region &&
           in_child_region == rhs.in_child_region;
}

bool pm::BlossomShatterEventData::operator!=(const pm::BlossomShatterEventData &rhs) const {
    return !(rhs == *this);
}

pm::MwpmEvent::MwpmEvent(pm::RegionHitRegionEventData region_hit_region_event_data)
    : region_hit_region_event_data(region_hit_region_event_data), event_type(pm::REGION_HIT_REGION) {
}
pm::MwpmEvent::MwpmEvent(pm::RegionHitBoundaryEventData region_hit_region_event_data)
    : region_hit_boundary_event_data(region_hit_region_event_data), event_type(pm::REGION_HIT_BOUNDARY) {
}
pm::MwpmEvent::MwpmEvent(pm::BlossomShatterEventData region_hit_region_event_data)
    : blossom_shatter_event_data(region_hit_region_event_data), event_type(pm::BLOSSOM_SHATTER) {
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
