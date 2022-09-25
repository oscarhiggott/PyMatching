#include "pymatching/fill_match/tracker/flood_check_event.h"

using namespace pm;

FloodCheckEvent::FloodCheckEvent(DetectorNode *data_look_at_node, cyclic_time_int time)
    : data_look_at_node(data_look_at_node), time(time), tentative_event_type(LOOK_AT_NODE) {
}

FloodCheckEvent::FloodCheckEvent(GraphFillRegion *data_look_at_shrinking_region, cyclic_time_int time)
    : data_look_at_shrinking_region(data_look_at_shrinking_region),
      time(time),
      tentative_event_type(LOOK_AT_SHRINKING_REGION) {
}

FloodCheckEvent::FloodCheckEvent(SearchDetectorNode *data_look_at_search_node, cyclic_time_int time)
        : data_look_at_search_node(data_look_at_search_node), time(time), tentative_event_type(LOOK_AT_SEARCH_NODE) {
}

FloodCheckEvent::FloodCheckEvent(cyclic_time_int time) : time(time), tentative_event_type(NO_FLOOD_CHECK_EVENT) {
}

bool FloodCheckEvent::operator==(const FloodCheckEvent &rhs) const {
    if (time != rhs.time || tentative_event_type != rhs.tentative_event_type)
        return false;
    switch (tentative_event_type) {
        case LOOK_AT_NODE:
            return data_look_at_node == rhs.data_look_at_node;
        case LOOK_AT_SHRINKING_REGION:
            return data_look_at_shrinking_region == rhs.data_look_at_shrinking_region;
        case LOOK_AT_SEARCH_NODE:
            return data_look_at_search_node == rhs.data_look_at_search_node;
        case NO_FLOOD_CHECK_EVENT:
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
            out << "LOOK_AT_SHRINKING_REGION, .region=" << ev.data_look_at_shrinking_region;
            break;
        case LOOK_AT_NODE:
            out << "LOOK_AT_NODE, .node=" << ev.data_look_at_node;
            break;
        case LOOK_AT_SEARCH_NODE:
            out << "LOOK_AT_SEARCH_NODE, .node=" << ev.data_look_at_search_node;
            break;
        case NO_FLOOD_CHECK_EVENT:
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
