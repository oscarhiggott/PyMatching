// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PYMATCHING2_MWPM_H
#define PYMATCHING2_MWPM_H

#include "pymatching/sparse_blossom/flooder/graph_flooder.h"
#include "pymatching/sparse_blossom/matcher/alternating_tree.h"
#include "pymatching/sparse_blossom/search/search_flooder.h"

namespace pm {

struct AltTreeNode;

struct MatchingResult {
    pm::obs_int obs_mask;
    total_weight_int weight;
    MatchingResult();

    bool operator==(const MatchingResult& rhs) const;

    bool operator!=(const MatchingResult& rhs) const;

    MatchingResult(obs_int obs_mask, total_weight_int weight);

    MatchingResult& operator+=(const MatchingResult& rhs);
    MatchingResult operator+(const MatchingResult& rhs) const;
};

struct Mwpm {
    GraphFlooder flooder;
    Arena<AltTreeNode> node_arena;
    SearchFlooder search_flooder;

    Mwpm();
    explicit Mwpm(GraphFlooder flooder);
    Mwpm(GraphFlooder flooder, SearchFlooder search_flooder);
    Mwpm(Mwpm&& other) noexcept;
    Mwpm& operator=(Mwpm&& other) noexcept;

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
    GraphFillRegion* pair_and_shatter_subblossoms_and_extract_matches(GraphFillRegion* region, MatchingResult& res);
    MatchingResult shatter_blossom_and_extract_matches(GraphFillRegion* region);

    GraphFillRegion* pair_and_shatter_subblossoms_and_extract_match_edges(
        GraphFillRegion* region, std::vector<CompressedEdge>& match_edges);
    void shatter_blossom_and_extract_match_edges(GraphFillRegion* region, std::vector<CompressedEdge>& match_edges);
    void extract_paths_from_match_edges(
        const std::vector<CompressedEdge>& match_edges, uint8_t* obs_begin_ptr, pm::total_weight_int& weight);

    void verify_invariants() const;

    void create_detection_event(DetectorNode* node);
    void reset();
};
}  // namespace pm

#endif  // PYMATCHING2_MWPM_H
