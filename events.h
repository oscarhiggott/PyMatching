#ifndef PYMATCHING2_EVENTS_H
#define PYMATCHING2_EVENTS_H

#include "varying.h"
#include "graph.h"
#include "graph_fill_region.h"

namespace pm{

    struct TentativeNeighborInteractionEventData {
        DetectorNode* detector_node_1;
        size_t node_1_neighbor_index;
        DetectorNode* detector_node_2;
        size_t node_2_neighbor_index;
        TentativeNeighborInteractionEventData(DetectorNode* detector_node_1,
                                              size_t node_1_neighbor_index,
                                              DetectorNode* detector_node_2,
                                              size_t node_2_neighbor_index);
        TentativeNeighborInteractionEventData() = default;
    };

    struct TentativeRegionShrinkEventData {
        GraphFillRegion* region;
        TentativeRegionShrinkEventData() = default;
        explicit TentativeRegionShrinkEventData(GraphFillRegion* region);
    };

    enum TentativeType : uint8_t {INTERACTION, SHRINKING};

    struct TentativeEvent{
        union {
            TentativeNeighborInteractionEventData neighbor_interaction_event_data;
            TentativeRegionShrinkEventData region_shrink_event_data;
        };
        pm::time_int time;
        TentativeType tentative_event_type;
        bool is_invalidated;

        TentativeEvent(DetectorNode* detector_node_1, size_t node_1_neighbor_index,
                       DetectorNode* detector_node_2, size_t node_2_neighbor_index,
                       time_int time);
        TentativeEvent(pm::GraphFillRegion* region, time_int time);
        TentativeEvent() = default;

        bool operator<(const TentativeEvent &rhs) const;

        bool operator>(const TentativeEvent &rhs) const;

        bool operator<=(const TentativeEvent &rhs) const;

        bool operator>=(const TentativeEvent &rhs) const;

        void invalidate();
    };

    struct RegionHitRegionEventData {
        GraphFillRegion* region1;
        GraphFillRegion* region2;
        CompressedEdge edge;
    };

    struct RegionHitBoundaryEventData {
        GraphFillRegion* region;
        CompressedEdge edge;
    };

    struct BlossomImplodeEventData {
        GraphFillRegion* blossom_region;
        GraphFillRegion* in_parent_region;
        GraphFillRegion* in_child_region;
    };

    enum MwpmEventType : uint8_t {
        REGION_HIT_REGION,
        REGION_HIT_BOUNDARY,
        BLOSSOM_IMPLODE
    };

    struct MwpmEvent {
        union {
            RegionHitRegionEventData region_hit_region_event_data;
            RegionHitBoundaryEventData region_hit_boundary_event_data;
            BlossomImplodeEventData blossom_implode_event_data;
        };
        MwpmEventType event_type;
    };

    inline bool pm::TentativeEvent::operator<(const pm::TentativeEvent &rhs) const {
        return time < rhs.time;
    }

    inline bool pm::TentativeEvent::operator>(const pm::TentativeEvent &rhs) const {
        return rhs < *this;
    }

    inline bool pm::TentativeEvent::operator<=(const pm::TentativeEvent &rhs) const {
        return !(rhs < *this);
    }

    inline bool pm::TentativeEvent::operator>=(const pm::TentativeEvent &rhs) const {
        return !(*this < rhs);
    }

}

#endif //PYMATCHING2_EVENTS_H
