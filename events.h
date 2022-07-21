#ifndef PYMATCHING2_EVENTS_H
#define PYMATCHING2_EVENTS_H

#include "varying.h"
#include "graph.h"
#include "graph_fill_region.h"




struct TentativeNeighborInteractionEventData {
    DetectorNode* detector_node_1;
    int schedule_index_1; // rename to neighbor index or something
    DetectorNode* detector_node_2;
    int schedule_index_2;
};

struct TentativeRegionShrinkEventData {
    GraphFillRegion* region;
};

// declare enum for tentative event types (inherit uint8_t)

// For all structs, biggest fields first

struct TentativeEvent{
    pm::time_int time;
    bool is_invalidated;
    //enum here
    union{
        TentativeNeighborInteractionEventData neighbor_interaction_event_data;
        TentativeRegionShrinkEventData tentative_region_shrink_event_data;
    };
};

// Same for MwpmEvents

struct MwpmEvent {};

struct BlossomImplodeEvent : MwpmEvent {
    GraphFillRegion* blossom_region;
    GraphFillRegion* in_parent_region;
    GraphFillRegion* in_child_region;
};

struct RegionHitRegionEvent : MwpmEvent {
    GraphFillRegion* region1;
    GraphFillRegion* region2;
    CompressedEdge edge;
};

struct RegionHitBoundaryEvent : MwpmEvent {
    GraphFillRegion* region;
    CompressedEdge edge;
};


#endif //PYMATCHING2_EVENTS_H
