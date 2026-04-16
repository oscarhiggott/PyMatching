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

// [NEW] Tests for K-MWM functionality.

#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"

#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/flooder/graph.h"

// Helper: build a d=5 repetition code matching graph.
//
//  [boundary] --w0-- 0 --w1-- 1 --w2-- 2 --w3-- 3 --w4-- [boundary]
//                    ^                             ^
//                observable 0                   (no obs)
//
// 4 detector nodes (0..3), 1 observable on the left boundary edge.
pm::Mwpm make_repetition_code_d5(
    pm::weight_int w0 = 100, pm::weight_int w1 = 100, pm::weight_int w2 = 100,
    pm::weight_int w3 = 100, pm::weight_int w4 = 100) {
    pm::MatchingGraph graph(4, 1);
    graph.add_boundary_edge(0, w0, {0});
    graph.add_edge(0, 1, w1, {});
    graph.add_edge(1, 2, w2, {});
    graph.add_edge(2, 3, w3, {});
    graph.add_boundary_edge(3, w4, {});
    return pm::Mwpm(pm::GraphFlooder(std::move(graph)));
}

TEST(KMwm, Milestone1_MatchEdgesConsistentWithDecode) {
    auto mwpm = make_repetition_code_d5();
    std::vector<uint64_t> detection_events = {0, 3};

    // Decode without match_edges (original behavior)
    auto ref = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false);
    mwpm.reset();

    // Decode with match_edges
    std::vector<pm::CompressedEdge> match_edges;
    auto result = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false, &match_edges);

    ASSERT_EQ(result.weight, ref.weight);
    ASSERT_EQ(result.obs_mask, ref.obs_mask);
    ASSERT_GE(match_edges.size(), 1u);

    pm::obs_int combined_obs = 0;
    for (auto& edge : match_edges)
        combined_obs ^= edge.obs_mask;
    ASSERT_EQ(combined_obs, result.obs_mask);
}

TEST(KMwm, Milestone1_EmptySyndrome) {
    auto mwpm = make_repetition_code_d5();
    std::vector<uint64_t> detection_events = {};
    std::vector<pm::CompressedEdge> match_edges;

    auto result = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false, &match_edges);

    ASSERT_EQ(result.weight, 0);
    ASSERT_EQ(result.obs_mask, 0u);
    ASSERT_EQ(match_edges.size(), 0u);
}

TEST(KMwm, Milestone1_BoundaryMatch) {
    auto mwpm = make_repetition_code_d5();
    std::vector<uint64_t> detection_events = {0};

    auto ref = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false);
    mwpm.reset();

    std::vector<pm::CompressedEdge> match_edges;
    auto result = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false, &match_edges);

    ASSERT_EQ(result.weight, ref.weight);
    ASSERT_EQ(result.obs_mask, ref.obs_mask);
    ASSERT_EQ(match_edges.size(), 1u);
    ASSERT_TRUE(match_edges[0].loc_from == nullptr || match_edges[0].loc_to == nullptr);
}

TEST(KMwm, Milestone1_TwoPairs) {
    // d=9 repetition code (8 detectors):
    //  [boundary] --100-- 0 --100-- 1 --100-- ... --100-- 7 --100-- [boundary]
    //
    // Detection events at {0,1,4,5} → two nearby pairs → should produce 2 match edges.
    pm::MatchingGraph graph(8, 1);
    graph.add_boundary_edge(0, 100, {0});
    for (int i = 0; i < 7; i++)
        graph.add_edge(i, i + 1, 100, {});
    graph.add_boundary_edge(7, 100, {});
    auto mwpm = pm::Mwpm(pm::GraphFlooder(std::move(graph)));

    std::vector<uint64_t> detection_events = {0, 1, 4, 5};

    auto ref = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false);
    mwpm.reset();

    std::vector<pm::CompressedEdge> match_edges;
    auto result = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false, &match_edges);

    ASSERT_EQ(result.weight, ref.weight);
    ASSERT_EQ(result.obs_mask, ref.obs_mask);
    ASSERT_EQ(match_edges.size(), 2u);
}

TEST(KMwm, Milestone1_AsymmetricWeights) {
    auto mwpm = make_repetition_code_d5(10, 100, 100, 100, 200);
    std::vector<uint64_t> detection_events = {0};

    auto ref = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false);
    mwpm.reset();

    std::vector<pm::CompressedEdge> match_edges;
    auto result = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false, &match_edges);

    ASSERT_EQ(result.weight, ref.weight);
    ASSERT_EQ(result.weight, 10);
    ASSERT_EQ(result.obs_mask, ref.obs_mask);
}

TEST(KMwm, Milestone1_GraphWithBlossom) {
    // Triangle graph (odd cycle) that triggers blossom formation during solving.
    //
    //     0 ---10--- 1
    //      \        /
    //      12     10
    //   (obs 0)  /
    //        \ /
    //         2 ---10--- 3 ---50--- [boundary]
    //
    // Detection events at all 4 nodes. The triangle 0-1-2 forms a blossom.
    pm::MatchingGraph graph(4, 1);
    graph.add_edge(0, 1, 10, {});
    graph.add_edge(1, 2, 10, {});
    graph.add_edge(0, 2, 12, {0});
    graph.add_edge(2, 3, 10, {});
    graph.add_boundary_edge(3, 50, {});
    auto mwpm = pm::Mwpm(pm::GraphFlooder(std::move(graph)));

    std::vector<uint64_t> detection_events = {0, 1, 2, 3};

    auto ref = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false);
    mwpm.reset();

    std::vector<pm::CompressedEdge> match_edges;
    auto result = pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events, false, &match_edges);

    ASSERT_EQ(result.weight, ref.weight);
    ASSERT_EQ(result.obs_mask, ref.obs_mask);
    ASSERT_GE(match_edges.size(), 2u);

    pm::obs_int combined_obs = 0;
    for (auto& edge : match_edges)
        combined_obs ^= edge.obs_mask;
    ASSERT_EQ(combined_obs, result.obs_mask);
}
