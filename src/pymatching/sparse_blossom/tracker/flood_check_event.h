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

#ifndef PYMATCHING2_FLOOD_CHECK_EVENT_H
#define PYMATCHING2_FLOOD_CHECK_EVENT_H

#include "pymatching/sparse_blossom/ints.h"

namespace pm {

struct DetectorNode;
struct GraphFillRegion;
struct SearchDetectorNode;

enum FloodCheckEventType : uint8_t {
    /// A placeholder value indicating there was no event.
    NO_FLOOD_CHECK_EVENT,

    /// Indicates that an event may be happening at a detector node. The event could be:
    /// - The node's region growing into an empty neighbor.
    /// - The node's region colliding with an adjacent boundary.
    /// - The node's region colliding with an adjacent region.
    LOOK_AT_NODE,

    /// Indicates that a region-level event might be happening. The event could be:
    /// - The region shrinking enough that a detector node needs to be removed from it.
    /// - The region being a blossom and shrinking to the point where it must shatter.
    /// - The region shrinking to point and causing a degenerate collision between its neighbors.
    LOOK_AT_SHRINKING_REGION,

    /// Indicates that an event may be happening at a SearchDetectorNode during a Dijkstra search.
    /// The event could be:
    /// - The node's exploratory region growing into an empty neighbor.
    /// - The node's exploratory region colliding with an adjacent boundary.
    /// - The node's exploratory region colliding with an adjacent exploratory region.
    LOOK_AT_SEARCH_NODE,
};

struct FloodCheckEvent {
    union {
        DetectorNode *data_look_at_node;
        GraphFillRegion *data_look_at_shrinking_region;
        SearchDetectorNode *data_look_at_search_node;
    };
    cyclic_time_int time;
    FloodCheckEventType tentative_event_type;

    FloodCheckEvent(DetectorNode *data, cyclic_time_int time);
    FloodCheckEvent(GraphFillRegion *data, cyclic_time_int time);
    FloodCheckEvent(SearchDetectorNode *data, cyclic_time_int time);
    explicit FloodCheckEvent(cyclic_time_int time);
    FloodCheckEvent() = delete;

    bool operator==(const FloodCheckEvent &rhs) const;
    bool operator!=(const FloodCheckEvent &rhs) const;

    std::string str() const;
};

inline FloodCheckEvent::FloodCheckEvent(DetectorNode *data_look_at_node, cyclic_time_int time)
    : data_look_at_node(data_look_at_node), time(time), tentative_event_type(LOOK_AT_NODE) {
}

inline FloodCheckEvent::FloodCheckEvent(GraphFillRegion *data_look_at_shrinking_region, cyclic_time_int time)
    : data_look_at_shrinking_region(data_look_at_shrinking_region),
      time(time),
      tentative_event_type(LOOK_AT_SHRINKING_REGION) {
}

inline FloodCheckEvent::FloodCheckEvent(SearchDetectorNode *data_look_at_search_node, cyclic_time_int time)
    : data_look_at_search_node(data_look_at_search_node), time(time), tentative_event_type(LOOK_AT_SEARCH_NODE) {
}

inline FloodCheckEvent::FloodCheckEvent(cyclic_time_int time) : time(time), tentative_event_type(NO_FLOOD_CHECK_EVENT) {
}

std::ostream &operator<<(std::ostream &out, const FloodCheckEvent &c);

}  // namespace pm

#endif  // PYMATCHING2_FLOOD_CHECK_EVENT_H
