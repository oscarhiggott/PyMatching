// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "pymatching/sparse_blossom/tracker/flood_check_event.h"

using namespace pm;

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
