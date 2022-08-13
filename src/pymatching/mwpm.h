#ifndef PYMATCHING2_MWPM_H
#define PYMATCHING2_MWPM_H

#include "pymatching/graph_flooder.h"

namespace pm {

struct MatchingResult {
    pm::obs_int obs_mask;
    pm::time_int weight;
    MatchingResult();

    bool operator==(const MatchingResult& rhs) const;

    bool operator!=(const MatchingResult& rhs) const;

    MatchingResult(obs_int obs_mask, time_int weight);

    MatchingResult& operator+=(const MatchingResult& rhs);
    MatchingResult operator+(const MatchingResult& rhs) const;
};

class Mwpm {
   public:
    GraphFlooder flooder;
    // detection_events can go, instead stored locally in decoding method
    std::vector<DetectorNode*> detection_events;
    void process_event(const MwpmEvent& event);
    explicit Mwpm(GraphFlooder& flooder);
    void handle_blossom_shattering(const BlossomShatterEventData& event);
    void shatter_descendants_into_matches_and_freeze(AltTreeNode& alt_tree_node);
    void handle_tree_hitting_boundary(const RegionHitBoundaryEventData& event);
    void handle_tree_hitting_match(
        GraphFillRegion* unmatched_region,
        GraphFillRegion* matched_region,
        const CompressedEdge& unmatched_to_matched_edge);
    void handle_tree_hitting_boundary_match(
        GraphFillRegion* unmatched_region,
        GraphFillRegion* matched_region,
        const CompressedEdge& unmatched_to_matched_edge);
    void handle_tree_hitting_self(const RegionHitRegionEventData& event, AltTreeNode* common_ancestor);
    void handle_tree_hitting_other_tree(const RegionHitRegionEventData& event);
    MatchingResult shatter_blossom_and_extract_matches(GraphFillRegion* region);
};
}  // namespace pm

#endif  // PYMATCHING2_MWPM_H
