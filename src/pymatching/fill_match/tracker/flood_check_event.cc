#include "pymatching/fill_match/tracker/flood_check_event.h"

using namespace pm;

FloodCheckEvent::FloodCheckEvent(TentativeEventData_LookAtNode data_look_at_node, cyclic_time_int time)
    : data_look_at_node(data_look_at_node),
      time(time),
      tentative_event_type(LOOK_AT_NODE) {
}

FloodCheckEvent::FloodCheckEvent(TentativeEventData_LookAtShrinkingRegion data_look_at_shrinking_region, cyclic_time_int time)
    : data_look_at_shrinking_region(data_look_at_shrinking_region), time(time), tentative_event_type(LOOK_AT_SHRINKING_REGION) {
}
FloodCheckEvent::FloodCheckEvent(cyclic_time_int time) : time(time), tentative_event_type(NO_TENTATIVE_EVENT) {
}


bool TentativeEventData_LookAtNode::operator==(const TentativeEventData_LookAtNode &rhs) const {
    return detector_node == rhs.detector_node;
}

bool TentativeEventData_LookAtNode::operator!=(const TentativeEventData_LookAtNode &rhs) const {
    return !(rhs == *this);
}

bool TentativeEventData_LookAtShrinkingRegion::operator==(const TentativeEventData_LookAtShrinkingRegion &rhs) const {
    return region == rhs.region;
}

bool TentativeEventData_LookAtShrinkingRegion::operator!=(const TentativeEventData_LookAtShrinkingRegion &rhs) const {
    return !(rhs == *this);
}

bool FloodCheckEvent::operator==(const FloodCheckEvent &rhs) const {
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

std::ostream &pm::operator<<(std::ostream &out, const FloodCheckEvent &ev) {
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

std::string FloodCheckEvent::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

bool FloodCheckEvent::operator!=(const FloodCheckEvent &rhs) const {
    return !(rhs == *this);
}
