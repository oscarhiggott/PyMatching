#ifndef PYMATCHING2_EVENTS_H
#define PYMATCHING2_EVENTS_H

#include "varying.h"
#include "graph.h"
#include "graph_fill_region.h"

namespace pm{

    struct TentativeNeighborInteractionEventData {
        DetectorNode* detector_node_1;

        bool operator==(const TentativeNeighborInteractionEventData &rhs) const;

        bool operator!=(const TentativeNeighborInteractionEventData &rhs) const;

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

        bool operator==(const TentativeRegionShrinkEventData &rhs) const;

        bool operator!=(const TentativeRegionShrinkEventData &rhs) const;

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

        bool operator==(const TentativeEvent &rhs) const;

        bool operator!=(const TentativeEvent &rhs) const;

        bool operator>=(const TentativeEvent &rhs) const;

        void invalidate();
    };

    struct RegionHitRegionEventData {
        GraphFillRegion* region1;
        GraphFillRegion* region2;
        CompressedEdge edge;
        RegionHitRegionEventData() = default;

        bool operator==(const RegionHitRegionEventData &rhs) const;

        bool operator!=(const RegionHitRegionEventData &rhs) const;

        RegionHitRegionEventData(GraphFillRegion* region1, GraphFillRegion* region2, CompressedEdge edge);
    };

    struct RegionHitBoundaryEventData {
        GraphFillRegion* region;
        CompressedEdge edge;

        bool operator==(const RegionHitBoundaryEventData &rhs) const;

        bool operator!=(const RegionHitBoundaryEventData &rhs) const;

        RegionHitBoundaryEventData() = default;
        RegionHitBoundaryEventData(GraphFillRegion* region, CompressedEdge edge);
    };

    struct BlossomShatterEventData {
        GraphFillRegion* blossom_region;
        GraphFillRegion* in_parent_region;
        GraphFillRegion* in_child_region;

        bool operator==(const BlossomShatterEventData &rhs) const;

        bool operator!=(const BlossomShatterEventData &rhs) const;

        BlossomShatterEventData() = default;
        BlossomShatterEventData(GraphFillRegion* blossom_region, GraphFillRegion* in_parent_region,
                                GraphFillRegion* in_child_region);
    };

    enum MwpmEventType : uint8_t {
        REGION_HIT_REGION,
        REGION_HIT_BOUNDARY,
        BLOSSOM_IMPLODE,
        NO_EVENT
    };

    struct MwpmEvent {
        union {
            RegionHitRegionEventData region_hit_region_event_data;
            RegionHitBoundaryEventData region_hit_boundary_event_data;
            BlossomShatterEventData blossom_shatter_event_data;
        };
        MwpmEventType event_type;
        MwpmEvent() = default;
        MwpmEvent(GraphFillRegion* region1, GraphFillRegion* region2, CompressedEdge edge);

        bool operator==(const MwpmEvent &rhs) const;

        bool operator!=(const MwpmEvent &rhs) const;

        MwpmEvent(GraphFillRegion* region, CompressedEdge edge);
        MwpmEvent(GraphFillRegion* blossom_region, GraphFillRegion* in_parent_region,
                  GraphFillRegion* in_child_region);
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
