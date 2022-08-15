#include "pymatching/events.h"

pm::TentativeEvent::TentativeEvent(pm::TentativeEventData_LookAtNode data_look_at_node, cyclic_time_int time)
    : data_look_at_node(data_look_at_node),
      time(time),
      tentative_event_type(LOOK_AT_NODE) {
}

pm::TentativeEvent::TentativeEvent(pm::TentativeEventData_LookAtShrinkingRegion data_look_at_shrinking_region, cyclic_time_int time)
    : data_look_at_shrinking_region(data_look_at_shrinking_region), time(time), tentative_event_type(LOOK_AT_SHRINKING_REGION) {
}
pm::TentativeEvent::TentativeEvent(cyclic_time_int time) : time(time), tentative_event_type(NO_TENTATIVE_EVENT) {
}


bool pm::TentativeEventData_LookAtNode::operator==(const pm::TentativeEventData_LookAtNode &rhs) const {
    return detector_node == rhs.detector_node;
}

bool pm::TentativeEventData_LookAtNode::operator!=(const pm::TentativeEventData_LookAtNode &rhs) const {
    return !(rhs == *this);
}

bool pm::TentativeEventData_LookAtShrinkingRegion::operator==(const pm::TentativeEventData_LookAtShrinkingRegion &rhs) const {
    return region == rhs.region;
}

bool pm::TentativeEventData_LookAtShrinkingRegion::operator!=(const pm::TentativeEventData_LookAtShrinkingRegion &rhs) const {
    return !(rhs == *this);
}

bool pm::TentativeEvent::operator==(const TentativeEvent &rhs) const {
    if (time != rhs.time || tentative_event_type != rhs.tentative_event_type)
        return false;
    switch (tentative_event_type) {
        case LOOK_AT_NODE:
            return data_look_at_node == rhs.data_look_at_node;
        case LOOK_AT_SHRINKING_REGION:
            return data_look_at_shrinking_region == rhs.data_look_at_shrinking_region;
        case NO_TENTATIVE_EVENT:
            return true;
        default:
            throw std::invalid_argument("Unrecognized event type");
    }
}

std::ostream &pm::operator<<(std::ostream &out, const TentativeEvent &ev) {
    out << "TentativeEvent{.time=";
    out << ev.time;
    out << ", .type=";
    switch (ev.tentative_event_type) {
        case LOOK_AT_SHRINKING_REGION:
            out << "LOOK_AT_SHRINKING_REGION, .region=" << ev.data_look_at_shrinking_region.region;
            break;
        case LOOK_AT_NODE:
            out << "LOOK_AT_NODE, .node=" << ev.data_look_at_node.detector_node;
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

std::ostream &pm::operator<<(std::ostream &out, const RegionHitBoundaryEventData &ev) {
    out << "{.region=" << ev.region;
    out << ", .edge=" << ev.edge;
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
            out << "REGION_HIT_BOUNDARY, .dat=" << ev.region_hit_boundary_event_data;
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

std::string pm::RegionHitBoundaryEventData::str() const {
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

bool pm::RegionHitRegionEventData::operator==(const pm::RegionHitRegionEventData &rhs) const {
    if (region1 == rhs.region1 && region2 == rhs.region2 && edge == rhs.edge) {
        return true;
    }

    // Also equal if reversed.
    return region2 == rhs.region1 && region1 == rhs.region2 && edge == rhs.edge.reversed();
}

bool pm::RegionHitRegionEventData::operator!=(const pm::RegionHitRegionEventData &rhs) const {
    return !(rhs == *this);
}

bool pm::RegionHitBoundaryEventData::operator==(const pm::RegionHitBoundaryEventData &rhs) const {
    return region == rhs.region && edge == rhs.edge;
}

bool pm::RegionHitBoundaryEventData::operator!=(const pm::RegionHitBoundaryEventData &rhs) const {
    return !(rhs == *this);
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
