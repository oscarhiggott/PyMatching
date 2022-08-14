#ifndef PYMATCHING2_EVENTS_H
#define PYMATCHING2_EVENTS_H

#include "pymatching/graph.h"
#include "pymatching/graph_fill_region.h"
#include "pymatching/varying.h"

namespace pm {

struct TentativeNeighborInteractionEventData {
    DetectorNode *detector_node_1;
    size_t node_1_neighbor_index;
    DetectorNode *detector_node_2;
    size_t node_2_neighbor_index;

    bool operator==(const TentativeNeighborInteractionEventData &rhs) const;
    bool operator!=(const TentativeNeighborInteractionEventData &rhs) const;
};

struct TentativeRegionShrinkEventData {
    GraphFillRegion *region;

    bool operator==(const TentativeRegionShrinkEventData &rhs) const;
    bool operator!=(const TentativeRegionShrinkEventData &rhs) const;
};

enum TentativeType : uint8_t {
    NO_TENTATIVE_EVENT,
    INTERACTION,
    SHRINKING
};

struct TentativeEvent {
    union {
        TentativeNeighborInteractionEventData neighbor_interaction_event_data;
        TentativeRegionShrinkEventData region_shrink_event_data;
    };
    pm::time_int time;
    TentativeType tentative_event_type;

    /// Validation index for the event. When events are created, this is set to a new unique value
    /// and the objects affected by the event are marked with the same unique value. When the event
    /// is actually processed, this value is compared to the markings on the affects objects. If
    /// they differ, the event is invalid.
    ///
    /// This makes invalidating events, starting from an object affected by that event, very cheap:
    /// simply increment its validation index. (Decrementing is not a safe way to invalidate because
    /// the previous index may have also been an event affecting that object.)
    uint64_t vid;

    TentativeEvent(
        TentativeNeighborInteractionEventData data,
        time_int time,
        uint64_t validation_index);
    TentativeEvent(TentativeRegionShrinkEventData data, time_int time, uint64_t validation_index);
    TentativeEvent(time_int time, uint64_t validation_index = 0);
    TentativeEvent() = default;

    bool operator==(const TentativeEvent &rhs) const;
    bool operator!=(const TentativeEvent &rhs) const;

    bool is_still_valid() const;

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
