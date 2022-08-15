#ifndef PYMATCHING2_EVENTS_H
#define PYMATCHING2_EVENTS_H

#include "pymatching/fill_match/flooder/varying.h"
#include "pymatching/fill_match/flooder/compressed_edge.h"

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
    /// One of the colliding regions.
    GraphFillRegion *region1;
    /// The other colliding region.
    GraphFillRegion *region2;
    /// The path between the two regions, anchored to detection events in the regions.
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
    /// The growing region that hit the boundary.
    GraphFillRegion *region;
    /// The path from the region to the boundary, anchored to a detection event in the region.
    CompressedEdge edge;

    bool operator==(const RegionHitBoundaryEventData &rhs) const;
    bool operator!=(const RegionHitBoundaryEventData &rhs) const;
    std::string str() const;
};
std::ostream &operator<<(std::ostream &out, const RegionHitBoundaryEventData &ev);

struct BlossomShatterEventData {
    /// The shrinking blossom region that has become empty and needs to be destroyed.
    GraphFillRegion *blossom_region;
    /// The blossom child that is linked to the alternating tree parent of the blossom.
    GraphFillRegion *in_parent_region;
    /// The blossom child that is linked to the alternating tree child of the blossom.
    GraphFillRegion *in_child_region;

    bool operator==(const BlossomShatterEventData &rhs) const;
    bool operator!=(const BlossomShatterEventData &rhs) const;
    std::string str() const;
};
std::ostream &operator<<(std::ostream &out, const BlossomShatterEventData &ev);

enum MwpmEventType : uint8_t { NO_EVENT, REGION_HIT_REGION, REGION_HIT_BOUNDARY, BLOSSOM_SHATTER };

/// A MwpmEvent is an interaction that the min-weight-perform-matching algorithm must react to.
///
/// An example of a tentative event that the flooder produces while running that DOES NOT produce
/// a MwpmEvent is a region growing enough to cover a new empty detector node.
///
/// An example of a tentative event tha the flooder produce while running that does produce a
/// MwpmEvent is a region growing into a detector node owned by another region, as this indicates
/// the regions have collided and a match (or blossom or etc) needs to be formed.
struct MwpmEvent {
    /// Event data. The `event_type` field determines which is populated.
    union {
        RegionHitRegionEventData region_hit_region_event_data;
        RegionHitBoundaryEventData region_hit_boundary_event_data;
        BlossomShatterEventData blossom_shatter_event_data;
    };
    /// Indicates the type of notification being sent to the mwpm algorithm.
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
