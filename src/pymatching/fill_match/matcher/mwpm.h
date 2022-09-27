#ifndef PYMATCHING2_MWPM_H
#define PYMATCHING2_MWPM_H

#include "pymatching/fill_match/flooder/graph_flooder.h"
#include "pymatching/fill_match/matcher/alternating_tree.h"
#include "pymatching/fill_match/search/search_flooder.h"

namespace pm {

struct AltTreeNode;

struct MatchingResult {
    pm::obs_int obs_mask;
    cumulative_time_int weight;
    MatchingResult();

    bool operator==(const MatchingResult& rhs) const;

    bool operator!=(const MatchingResult& rhs) const;

    MatchingResult(obs_int obs_mask, cumulative_time_int weight);

    MatchingResult& operator+=(const MatchingResult& rhs);
    MatchingResult operator+(const MatchingResult& rhs) const;
};


struct ExtendedMatchingResult {
    std::vector<uint8_t> obs_crossed;
    cumulative_time_int weight;
    ExtendedMatchingResult();
    explicit ExtendedMatchingResult(size_t num_observables);

    bool operator==(const ExtendedMatchingResult& rhs) const;

    bool operator!=(const ExtendedMatchingResult& rhs) const;

    ExtendedMatchingResult(std::vector<uint8_t> obs_crossed, cumulative_time_int weight);

    ExtendedMatchingResult& operator+=(const ExtendedMatchingResult& rhs);
    ExtendedMatchingResult operator+(const ExtendedMatchingResult& rhs) const;
};


struct Mwpm {
    GraphFlooder flooder;
    Arena<AltTreeNode> node_arena;
    SearchFlooder search_flooder;

    explicit Mwpm(GraphFlooder flooder);
    Mwpm(GraphFlooder flooder, SearchFlooder search_flooder);

    AltTreeNode* make_child(
        AltTreeNode& parent,
        GraphFillRegion* child_inner_region,
        GraphFillRegion* child_outer_region,
        const CompressedEdge& child_inner_to_outer_edge,
        const CompressedEdge& child_compressed_edge);
    void process_event(const pm::MwpmEvent& event);
    void handle_blossom_shattering(const BlossomShatterEventData& event);
    void shatter_descendants_into_matches_and_freeze(AltTreeNode& alt_tree_node);
    void handle_tree_hitting_boundary(const RegionHitBoundaryEventData& event);
    void handle_region_hit_region(const MwpmEvent event);
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
    GraphFillRegion* pair_and_shatter_subblossoms_and_extract_matches(GraphFillRegion* region, MatchingResult &res);
    MatchingResult shatter_blossom_and_extract_matches(GraphFillRegion* region);

    GraphFillRegion* pair_and_shatter_subblossoms_and_extract_match_edges(
            GraphFillRegion* region,
            std::vector<CompressedEdge>& match_edges
            );
    void shatter_blossom_and_extract_match_edges(GraphFillRegion* region, std::vector<CompressedEdge>& match_edges);
    void extract_paths_from_match_edges(
            const std::vector<CompressedEdge>& match_edges,
            uint8_t *obs_begin_ptr,
            pm::cumulative_time_int& weight
            );

    void verify_invariants() const;

    void create_detection_event(DetectorNode* node);
};
}  // namespace pm

#endif  // PYMATCHING2_MWPM_H
