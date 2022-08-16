#ifndef PYMATCHING2_FLOOD_CHECK_EVENT_H
#define PYMATCHING2_FLOOD_CHECK_EVENT_H

#include "pymatching/fill_match/ints.h"

namespace pm {

struct DetectorNode;
struct GraphFillRegion;

enum FloodCheckEventType : uint8_t {
    /// A placeholder value indicating there was no event.
    NO_FLOOD_CHECK_EVENT,

    /// Indicates that an event may be happening at a detector node. The event could be:
    /// - The node's region colliding with an adjacent boundary.
    /// - The node's region colliding with an adjacent region.
    LOOK_AT_NODE,

    /// Indicates that a region-level event might be happening. The event could be:
    /// - The region shrinking enough that a detector node needs to be removed from it.
    /// - The region being a blossom and shrinking to the point where it must shatter.
    /// - The region shrinking to point and causing a degenerate collision between its neighbors.
    LOOK_AT_SHRINKING_REGION,
};

struct FloodCheckEvent {
    union {
        DetectorNode *data_look_at_node;
        GraphFillRegion *data_look_at_shrinking_region;
    };
    cyclic_time_int time;
    FloodCheckEventType tentative_event_type;

    FloodCheckEvent(DetectorNode *data, cyclic_time_int time);
    FloodCheckEvent(GraphFillRegion *data, cyclic_time_int time);
    explicit FloodCheckEvent(cyclic_time_int time);
    FloodCheckEvent() = delete;

    bool operator==(const FloodCheckEvent &rhs) const;
    bool operator!=(const FloodCheckEvent &rhs) const;

    std::string str() const;
};

std::ostream &operator<<(std::ostream &out, const FloodCheckEvent &c);

}  // namespace pm

#endif  // PYMATCHING2_FLOOD_CHECK_EVENT_H
