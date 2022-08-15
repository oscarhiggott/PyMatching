#ifndef PYMATCHING2_EVENTS_H
#define PYMATCHING2_EVENTS_H

#include "pymatching/varying.h"
#include "pymatching/compressed_edge.h"

namespace pm {

struct DetectorNode;
struct GraphFillRegion;

enum TentativeType : uint8_t {
    /// A placeholder value indicating there was no event.
    NO_TENTATIVE_EVENT,

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

struct TentativeEventData_LookAtNode {
    DetectorNode *detector_node;
    bool operator==(const TentativeEventData_LookAtNode &rhs) const;
    bool operator!=(const TentativeEventData_LookAtNode &rhs) const;
};

struct TentativeEventData_LookAtShrinkingRegion {
    GraphFillRegion *region;
    bool operator==(const TentativeEventData_LookAtShrinkingRegion &rhs) const;
    bool operator!=(const TentativeEventData_LookAtShrinkingRegion &rhs) const;
};

struct TentativeEvent {
    union {
        TentativeEventData_LookAtNode data_look_at_node;
        TentativeEventData_LookAtShrinkingRegion data_look_at_shrinking_region;
    };
    cyclic_time_int time;
    TentativeType tentative_event_type;

    TentativeEvent(TentativeEventData_LookAtNode data, cyclic_time_int time);
    TentativeEvent(TentativeEventData_LookAtShrinkingRegion data, cyclic_time_int time);
    explicit TentativeEvent(cyclic_time_int time);
    TentativeEvent() = default;

    bool operator==(const TentativeEvent &rhs) const;
    bool operator!=(const TentativeEvent &rhs) const;

    std::string str() const;
};

std::ostream &operator<<(std::ostream &out, const TentativeEvent &c);


struct RegionHitRegionEventData {
    GraphFillRegion *region1;
    GraphFillRegion *region2;
    CompressedEdge edge;

    inline RegionHitRegionEventData reversed() const {
        return {region2, region1, edge.reversed()};
    }
    bool operator==(const RegionHitRegionEventData &rhs) const;
    bool operator!=(const RegionHitRegionEventData &rhs) const;
    std::string str() const;
};
std::ostream &operator<<(std::ostream &out, const RegionHitRegionEventData &dat);

struct RegionHitBoundaryEventData {
    GraphFillRegion *region;
    CompressedEdge edge;

    bool operator==(const RegionHitBoundaryEventData &rhs) const;
    bool operator!=(const RegionHitBoundaryEventData &rhs) const;
};
std::ostream &operator<<(std::ostream &out, const RegionHitBoundaryEventData &ev);

struct BlossomShatterEventData {
    GraphFillRegion *blossom_region;
    GraphFillRegion *in_parent_region;
    GraphFillRegion *in_child_region;

    bool operator==(const BlossomShatterEventData &rhs) const;
    bool operator!=(const BlossomShatterEventData &rhs) const;
    std::string str() const;
};
std::ostream &operator<<(std::ostream &out, const BlossomShatterEventData &ev);

enum MwpmEventType : uint8_t { NO_EVENT, REGION_HIT_REGION, REGION_HIT_BOUNDARY, BLOSSOM_SHATTER };

struct MwpmEvent {
    union {
        RegionHitRegionEventData region_hit_region_event_data;
        RegionHitBoundaryEventData region_hit_boundary_event_data;
        BlossomShatterEventData blossom_shatter_event_data;
    };
    MwpmEventType event_type;

    MwpmEvent();
    MwpmEvent(RegionHitRegionEventData region_hit_region_event_data);      // NOLINT(google-explicit-constructor)
    MwpmEvent(RegionHitBoundaryEventData region_hit_boundary_event_data);  // NOLINT(google-explicit-constructor)
    MwpmEvent(BlossomShatterEventData blossom_shatter_event_data);         // NOLINT(google-explicit-constructor)
    inline static MwpmEvent no_event() {
        return {};
    }

    bool operator==(const MwpmEvent &rhs) const;
    bool operator!=(const MwpmEvent &rhs) const;
    std::string str() const;
};
std::ostream &operator<<(std::ostream &out, const MwpmEvent &ev);

}  // namespace pm

#endif  // PYMATCHING2_EVENTS_H
