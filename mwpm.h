#ifndef PYMATCHING2_MWPM_H
#define PYMATCHING2_MWPM_H

#include "graph_flooder.h"


class MatchingResult {
    std::vector<CompressedEdge> match_edges;
    weight_int total_weight;
    obs_int observables;
};


class Mwpm {
public:
    GraphFlooder flooder;
    // detection_events can go, instead stored locally in decoding method
    std::vector<DetectorNode*> detection_events;
    void add_detection_event(int detector_node_id);
    void process_event(const MwpmEvent& event);
    Mwpm(GraphFlooder flooder);
    void handle_blossom_imploding(const BlossomImplodeEventData& event);
    void handle_tree_hitting_match(const RegionHitRegionEventData& event);
    void shatter_descendants_into_matches_and_freeze(AltTreeNode& alt_tree_node);
    void handle_tree_hitting_boundary(const RegionHitBoundaryEventData& event);
    void handle_tree_hitting_boundary_match(const RegionHitRegionEventData& event);
    void handle_tree_hitting_self(const RegionHitRegionEventData& event);
    void handle_tree_hitting_other_tree(const RegionHitRegionEventData& event);
    MatchingResult extract_matching_and_reset_graph();
};

#endif //PYMATCHING2_MWPM_H
