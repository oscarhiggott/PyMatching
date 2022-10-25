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

#include "pymatching/sparse_blossom/matcher/mwpm.h"

#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/arena.h"
#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"
#include "pymatching/sparse_blossom/matcher/alternating_tree.h"

using namespace pm;

struct MwpmAltTreeTestData {
    Mwpm& mwpm;
    explicit MwpmAltTreeTestData(Mwpm& mwpm) : mwpm(mwpm) {
    }
    AltTreeEdge t(
        std::vector<AltTreeEdge> children, size_t inner_region_id, size_t outer_region_id, bool root = false) {
        AltTreeNode* node = mwpm.node_arena.alloc_unconstructed();
        auto& ns = mwpm.flooder.graph.nodes;
        CompressedEdge parent_ce{nullptr, nullptr, 0};
        if (root) {
            new (node) AltTreeNode(ns[outer_region_id].region_that_arrived);
        } else {
            new (node) AltTreeNode(
                ns[inner_region_id].region_that_arrived,
                ns[outer_region_id].region_that_arrived,
                {&ns[inner_region_id], &ns[outer_region_id], 0});
        }
        auto edge = AltTreeEdge(node, parent_ce);
        for (auto& child : children) {
            child.edge.loc_from = &ns[outer_region_id];
            child.edge.loc_to = child.alt_tree_node->inner_to_outer_edge.loc_from;
            edge.alt_tree_node->add_child(child);
        }
        return edge;
    }
};

bool regions_matched(std::vector<DetectorNode>& ns, size_t i, size_t j) {
    bool regions_correct = ns[i].region_that_arrived->match.region == ns[j].region_that_arrived &&
                           ns[j].region_that_arrived->match.region == ns[i].region_that_arrived;
    bool edge_locs_correct = ns[i].region_that_arrived->match.edge.loc_from == &ns[i] &&
                             ns[i].region_that_arrived->match.edge.loc_to == &ns[j] &&
                             ns[j].region_that_arrived->match.edge.loc_from == &ns[j] &&
                             ns[j].region_that_arrived->match.edge.loc_to == &ns[i];
    return regions_correct && edge_locs_correct;
}

MwpmEvent rhr(std::vector<DetectorNode>& ns, size_t i, size_t j) {
    return RegionHitRegionEventData{ns[i].region_that_arrived, ns[j].region_that_arrived, {&ns[i], &ns[j], 0}};
}

MwpmEvent rhr(std::vector<DetectorNode>& ns, size_t i, size_t j, obs_int obs_mask) {
    return RegionHitRegionEventData{ns[i].region_that_arrived, ns[j].region_that_arrived, {&ns[i], &ns[j], obs_mask}};
}

TEST(Mwpm, BlossomCreatedThenShattered) {
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(10, 64)));
    auto& g = mwpm.flooder.graph;
    g.add_edge(0, 1, 10, {0});
    g.add_edge(1, 4, 20, {1});
    g.add_edge(4, 3, 20, {0, 1});
    g.add_edge(3, 2, 12, {2});
    g.add_edge(0, 2, 16, {0, 2});
    g.add_edge(4, 5, 50, {1, 2});
    g.add_edge(2, 6, 100, {0, 1, 2});
    g.add_boundary_edge(5, 36, {3});
    for (size_t i = 0; i < 7; i++) {
        mwpm.create_detection_event(&mwpm.flooder.graph.nodes[i]);
    }
    auto& ns = mwpm.flooder.graph.nodes;

    auto dequeued_and_processed = [&](size_t expected_time) {
        auto ev = mwpm.flooder.run_until_next_mwpm_notification();
        EXPECT_EQ(expected_time, mwpm.flooder.queue.cur_time);
        mwpm.process_event(ev);
        return ev;
    };

    // 0:1 and 2:3 pair up.
    ASSERT_EQ(
        dequeued_and_processed(5),
        MwpmEvent(RegionHitRegionEventData{ns[0].region_that_arrived, ns[1].region_that_arrived, {&ns[0], &ns[1], 1}}));
    ASSERT_EQ(
        dequeued_and_processed(6),
        MwpmEvent(RegionHitRegionEventData{ns[2].region_that_arrived, ns[3].region_that_arrived, {&ns[2], &ns[3], 4}}));

    // 4 hits frozen 3 creating 4-3-2 tree.
    ASSERT_EQ(
        dequeued_and_processed(14),
        MwpmEvent(RegionHitRegionEventData{ns[3].region_that_arrived, ns[4].region_that_arrived, {&ns[3], &ns[4], 3}}));
    // 4 hits frozen 1 creating 0-1-4-3-2 tree.
    ASSERT_EQ(
        dequeued_and_processed(15),
        MwpmEvent(RegionHitRegionEventData{ns[1].region_that_arrived, ns[4].region_that_arrived, {&ns[1], &ns[4], 2}}));

    // 0 hits 2 creating blossom (0-1-4-3-2).
    ASSERT_EQ(
        dequeued_and_processed(17),
        MwpmEvent(RegionHitRegionEventData{ns[0].region_that_arrived, ns[2].region_that_arrived, {&ns[0], &ns[2], 5}}));
    auto* blossom = ns[0].region_that_arrived->blossom_parent;
    ASSERT_NE(blossom, nullptr);
    ASSERT_EQ(blossom->blossom_children.size(), 5);
    std::vector<RegionEdge> expected_blossom_children = {
        {ns[2].region_that_arrived, {&ns[2], &ns[3], 4}},
        {ns[3].region_that_arrived, {&ns[3], &ns[4], 3}},
        {ns[4].region_that_arrived, {&ns[4], &ns[1], 2}},
        {ns[1].region_that_arrived, {&ns[1], &ns[0], 1}},
        {ns[0].region_that_arrived, {&ns[0], &ns[2], 5}},
    };
    ASSERT_EQ(blossom->blossom_children, expected_blossom_children);
    for (auto& c : blossom->blossom_children) {
        ASSERT_TRUE(c.region->alt_tree_node == nullptr);
    }

    // 4 hits 5 creating matches
    ASSERT_EQ(
        dequeued_and_processed(25),
        MwpmEvent(RegionHitRegionEventData{blossom, ns[5].region_that_arrived, {&ns[4], &ns[5], 6}}));
    ASSERT_EQ(blossom->radius, VaryingCT(8 << 2));

    // Region 6 collides with matched blossom
    ASSERT_EQ(
        dequeued_and_processed(83),
        MwpmEvent(RegionHitRegionEventData{blossom, ns[6].region_that_arrived, {&ns[2], &ns[6], 7}}));
    ASSERT_EQ(ns[6].region_that_arrived->alt_tree_node->children[0].alt_tree_node->inner_region, blossom);

    // Blossom shatters
    ASSERT_EQ(
        dequeued_and_processed(91),
        MwpmEvent(BlossomShatterEventData{blossom, ns[2].region_that_arrived, ns[4].region_that_arrived}));
    ASSERT_EQ(ns[0].region_that_arrived->match.region, ns[1].region_that_arrived);
    ASSERT_EQ(ns[1].region_that_arrived->match.region, ns[0].region_that_arrived);

    auto e9 = mwpm.flooder.run_until_next_mwpm_notification();
    ASSERT_EQ(mwpm.flooder.queue.cur_time, 94);
    ASSERT_EQ(e9.event_type, REGION_HIT_BOUNDARY);
    mwpm.process_event(e9);
    auto e10 = mwpm.flooder.run_until_next_mwpm_notification();
    ASSERT_EQ(e10.event_type, NO_EVENT);
    ASSERT_EQ(ns[3].region_that_arrived->blossom_parent, nullptr);
    ASSERT_TRUE(regions_matched(ns, 4, 3));
    ASSERT_TRUE(regions_matched(ns, 2, 6));
    ASSERT_TRUE(regions_matched(ns, 1, 0));
    ASSERT_EQ(ns[5].region_that_arrived->match, (Match{nullptr, {&ns[5], nullptr, 8}}));
}

TEST(Mwpm, BlossomShatterDrivenWithoutFlooder) {
    size_t n = 10;
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(n + 3, 64)));
    auto& ns = mwpm.flooder.graph.nodes;
    for (size_t i = 0; i < n + 3; i++) {
        mwpm.create_detection_event(&ns[i]);
    }

    // Pair up
    for (size_t i = 0; i < n; i += 2) {
        mwpm.process_event(rhr(ns, i, i + 1));
    }
    // Form an alternating path
    for (size_t i = 1; i < n; i += 2) {
        mwpm.process_event(rhr(ns, n - i, n - i + 1));
    }
    // Close the path into a blossom
    mwpm.process_event(rhr(ns, 0, n));
    auto blossom = ns[0].region_that_arrived->blossom_parent;
    // Make the blossom become an inner node
    mwpm.process_event(RegionHitRegionEventData{ns[n + 1].region_that_arrived, blossom, {&ns[n + 1], &ns[5], 0}});
    mwpm.process_event(RegionHitRegionEventData{ns[n + 2].region_that_arrived, blossom, {&ns[n + 2], &ns[10], 0}});
    ASSERT_EQ(blossom->blossom_children.size(), 11);
    mwpm.process_event(BlossomShatterEventData{blossom, ns[10].region_that_arrived, ns[5].region_that_arrived});
    ASSERT_TRUE(regions_matched(ns, 8, 9));
    ASSERT_TRUE(regions_matched(ns, 6, 7));
    auto actual_tree = ns[12].region_that_arrived->alt_tree_node;
    ASSERT_EQ(actual_tree->parent.alt_tree_node, nullptr);

    MwpmAltTreeTestData d(mwpm);
    auto expected_tree = d.t({d.t({d.t({d.t({d.t({}, 5, 11)}, 3, 4)}, 1, 2)}, 10, 0)}, -1, 12, true).alt_tree_node;
    ASSERT_EQ(*actual_tree, *expected_tree);
}

TEST(Mwpm, BranchingTreeFormsBlossomThenHitsBoundaryMatch) {
    size_t n = 14;
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(n, 64)));
    auto& ns = mwpm.flooder.graph.nodes;
    for (size_t i = 0; i < n; i++) {
        mwpm.create_detection_event(&ns[i]);
    }

    // Form matches
    mwpm.process_event(rhr(ns, 1, 3));
    mwpm.process_event(rhr(ns, 2, 4));
    mwpm.process_event(rhr(ns, 5, 9));
    mwpm.process_event(rhr(ns, 6, 10));
    mwpm.process_event(rhr(ns, 7, 11));
    mwpm.process_event(rhr(ns, 8, 12));
    ASSERT_TRUE(regions_matched(ns, 1, 3));
    ASSERT_TRUE(regions_matched(ns, 2, 4));
    ASSERT_TRUE(regions_matched(ns, 5, 9));
    ASSERT_TRUE(regions_matched(ns, 6, 10));
    ASSERT_TRUE(regions_matched(ns, 7, 11));
    ASSERT_TRUE(regions_matched(ns, 8, 12));
    // Form alternating tree
    mwpm.process_event(rhr(ns, 0, 1));
    mwpm.process_event(rhr(ns, 0, 2));
    mwpm.process_event(rhr(ns, 3, 5));
    mwpm.process_event(rhr(ns, 3, 6));
    mwpm.process_event(rhr(ns, 4, 7));
    mwpm.process_event(rhr(ns, 4, 8));
    MwpmAltTreeTestData d(mwpm);
    auto expected_tree =
        d.t({d.t({d.t({}, 5, 9), d.t({}, 6, 10)}, 1, 3), d.t({d.t({}, 7, 11), d.t({}, 8, 12)}, 2, 4)}, -1, 0, true)
            .alt_tree_node;
    auto actual_tree = ns[0].region_that_arrived->alt_tree_node;
    ASSERT_EQ(*expected_tree, *actual_tree);
    std::vector<AltTreeEdge> expected_blossom_tree_children = {
        expected_tree->children[0].alt_tree_node->children[0], expected_tree->children[1].alt_tree_node->children[1]};
    // Form blossom
    mwpm.process_event(rhr(ns, 10, 11));
    auto blossom = ns[0].region_that_arrived->blossom_parent;
    ASSERT_EQ(blossom->blossom_children.size(), 9);
    ASSERT_EQ(blossom->alt_tree_node->children, expected_blossom_tree_children);
    // Region 13 matches to boundary
    mwpm.process_event(RegionHitBoundaryEventData{ns[13].region_that_arrived, {&ns[13], nullptr, 0}});
    // Blossom matches to region 13
    mwpm.process_event(MwpmEvent(RegionHitRegionEventData{blossom, ns[13].region_that_arrived, {&ns[4], &ns[13], 1}}));
    ASSERT_TRUE(regions_matched(ns, 5, 9));
    ASSERT_TRUE(regions_matched(ns, 8, 12));
    ASSERT_EQ(blossom->match, (Match{ns[13].region_that_arrived, {&ns[4], &ns[13], 1}}));
    ASSERT_EQ(ns[13].region_that_arrived->match, (Match{blossom, {&ns[13], &ns[4], 1}}));
}

TEST(Mwpm, BoundaryMatchHitsTree) {
    size_t n = 6;
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(n, 64)));
    auto& ns = mwpm.flooder.graph.nodes;
    for (size_t i = 0; i < n; i++)
        mwpm.create_detection_event(&ns[i]);
    // Form alternating tree
    mwpm.process_event(rhr(ns, 1, 2));
    mwpm.process_event(rhr(ns, 3, 4));
    mwpm.process_event(rhr(ns, 0, 1));
    mwpm.process_event(rhr(ns, 0, 3));
    // 5 matches to boundary
    mwpm.process_event(RegionHitBoundaryEventData{ns[5].region_that_arrived, {&ns[5], nullptr, 1}});
    // Tree matches to 5
    mwpm.process_event(rhr(ns, 5, 2));
    ASSERT_TRUE(regions_matched(ns, 5, 2));
    ASSERT_TRUE(regions_matched(ns, 1, 0));
    ASSERT_TRUE(regions_matched(ns, 3, 4));
}

std::vector<size_t> obs_mask_to_set_bits(pm::obs_int obs_mask) {
    std::vector<size_t> out;
    for (size_t i = 0; i < sizeof(pm::obs_int) * 8; i++) {
        if (obs_mask & ((pm::obs_int)1 << i))
            out.push_back(i);
    }
    return out;
}

TEST(Mwpm, ShatterBlossomAndExtractMatchesForPair) {
    size_t num_nodes = 20;
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(num_nodes, 64)));
    auto& g = mwpm.flooder.graph;
    for (size_t i = 0; i < num_nodes - 1; i++)
        g.add_edge(i, i + 1, 2, obs_mask_to_set_bits(i));
    auto& ns = mwpm.flooder.graph.nodes;
    mwpm.create_detection_event(&ns[6]);
    mwpm.create_detection_event(&ns[13]);
    auto e = mwpm.flooder.run_until_next_mwpm_notification();
    mwpm.process_event(e);
    ASSERT_EQ(ns[6].region_that_arrived->match.region, ns[13].region_that_arrived);
    auto res = mwpm.shatter_blossom_and_extract_matches(ns[6].region_that_arrived);
    ASSERT_EQ(res, MatchingResult(6 ^ 7 ^ 8 ^ 9 ^ 10 ^ 11 ^ 12, 14));
    for (auto& n : ns) {
        ASSERT_EQ(n.region_that_arrived, nullptr);
        ASSERT_EQ(n.radius_of_arrival, 0);
        ASSERT_EQ(n.observables_crossed_from_source, 0);
        ASSERT_EQ(n.reached_from_source, nullptr);
    }
}

obs_int set_bits_to_obs_mask(const std::vector<size_t>& set_bits) {
    obs_int obs_mask = 0;
    for (auto i : set_bits)
        obs_mask ^= (1 << i);
    return obs_mask;
}

TEST(Mwpm, ObsMaskToSetBits) {
    for (size_t i = 0; i < 64; i++) {
        auto x = obs_mask_to_set_bits(i);
        auto y = set_bits_to_obs_mask(x);
        ASSERT_EQ(y, i);
    }
}

void form_blossom_and_nested_blossom_then_match_example(Mwpm& mwpm) {
    auto& ns = mwpm.flooder.graph.nodes;
    cumulative_time_int w = 0;
    for (auto& n : ns) {
        w += 1;
        mwpm.create_detection_event(&n);
        n.region_that_arrived->radius = Varying((w << 2) + 1);
    }
    // Form first blossom
    mwpm.process_event(rhr(ns, 0, 1, 1 << 1));
    mwpm.process_event(rhr(ns, 2, 3, 1 << 2));
    mwpm.process_event(rhr(ns, 4, 3, 1 << 3));
    mwpm.process_event(rhr(ns, 2, 1, 1 << 4));
    mwpm.process_event(rhr(ns, 0, 4, 1 << 5));
    auto blossom1 = ns[0].region_that_arrived->blossom_parent;
    blossom1->radius = blossom1->radius + 2;
    ASSERT_EQ(blossom1->blossom_children.size(), 5);
    // Form second blossom
    mwpm.process_event(rhr(ns, 6, 7, 1 << 6));
    mwpm.process_event(rhr(ns, 5, 6, 1 << 7));
    mwpm.process_event(rhr(ns, 5, 7, 1 << 8));
    auto blossom2 = ns[6].region_that_arrived->blossom_parent;
    blossom2->radius = blossom2->radius + 3;
    ASSERT_EQ(blossom2->blossom_children.size(), 3);
    // Form third blossom
    mwpm.process_event(rhr(ns, 10, 11, 1 << 9));
    mwpm.process_event(rhr(ns, 8, 10, 1 << 10));
    mwpm.process_event(rhr(ns, 8, 11, 1 << 11));
    auto blossom3 = ns[10].region_that_arrived->blossom_parent;
    blossom3->radius = blossom3->radius + 5;
    ASSERT_EQ(blossom3->blossom_children.size(), 3);
    // Fourth blossom, containing blossom2, blossom3 and region 9
    mwpm.process_event(RegionHitRegionEventData{blossom2, blossom3, {&ns[5], &ns[8], 1 << 12}});
    mwpm.process_event(RegionHitRegionEventData{ns[9].region_that_arrived, blossom2, {&ns[9], &ns[5], 1 << 13}});
    mwpm.process_event(RegionHitRegionEventData{ns[9].region_that_arrived, blossom3, {&ns[9], &ns[8], 1 << 14}});
    auto blossom4 = blossom3->blossom_parent;
    blossom4->radius = blossom4->radius + 11;
    ASSERT_EQ(blossom4->blossom_children.size(), 3);
    // blossom1 matches to blossom2
    mwpm.process_event(RegionHitRegionEventData{blossom1, blossom4, {&ns[0], &ns[5], 1 << 15}});
}

TEST(Mwpm, ShatterBlossomAndExtractMatches) {
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(12, 64)));
    form_blossom_and_nested_blossom_then_match_example(mwpm);
    auto& ns = mwpm.flooder.graph.nodes;

    auto res = mwpm.shatter_blossom_and_extract_matches(ns[0].region_that_arrived->blossom_parent);

    ASSERT_EQ(res.obs_mask, set_bits_to_obs_mask({3, 4, 15, 6, 14, 9}));
    ASSERT_EQ(res.weight, 2 + 3 + 4 + 5 + 1 + 2 + 6 + 3 + 11 + 7 + 8 + 9 + 10 + 5 + 11 + 12);
}

TEST(Mwpm, ShatterBlossomAndExtractMatchEdges) {
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(12, 64)));
    form_blossom_and_nested_blossom_then_match_example(mwpm);
    auto& ns = mwpm.flooder.graph.nodes;

    std::vector<CompressedEdge> match_edges;
    mwpm.shatter_blossom_and_extract_match_edges(ns[0].region_that_arrived->blossom_parent, match_edges);

    pm::obs_int obs_mask = 0;
    for (auto& e : match_edges) {
        obs_mask ^= e.obs_mask;
    }

    ASSERT_EQ(obs_mask, set_bits_to_obs_mask({3, 4, 15, 6, 14, 9}));
    ASSERT_EQ(match_edges.size(), 6);

    std::vector<CompressedEdge> expected_edges = {
        {&ns[4], &ns[3], 1 << 3},
        {&ns[2], &ns[1], 1 << 4},
        {&ns[11], &ns[10], 1 << 9},
        {&ns[9], &ns[8], 1 << 14},
        {&ns[7], &ns[6], 1 << 6},
        {&ns[0], &ns[5], 1 << 15}};
    ASSERT_EQ(expected_edges, match_edges);
}

void form_nested_blossom_and_match_to_boundary(Mwpm& mwpm) {
    auto& ns = mwpm.flooder.graph.nodes;
    cumulative_time_int w = 0;
    for (auto& n : ns) {
        w += 1;
        mwpm.create_detection_event(&n);
        n.region_that_arrived->radius = Varying((w << 2) + 1);
    }
    // Form first blossom
    mwpm.process_event(rhr(ns, 1, 2, 1 << 1));
    mwpm.process_event(rhr(ns, 0, 1, 1 << 2));
    mwpm.process_event(rhr(ns, 0, 2, 1 << 3));
    auto blossom1 = ns[0].region_that_arrived->blossom_parent;
    ASSERT_EQ(blossom1->blossom_children.size(), 3);
    blossom1->radius = blossom1->radius + 20;
    mwpm.process_event(rhr(ns, 3, 4, 1 << 4));
    mwpm.process_event(RegionHitRegionEventData{blossom1, ns[3].region_that_arrived, {&ns[0], &ns[3], 1 << 10}});
    mwpm.process_event(RegionHitRegionEventData{blossom1, ns[4].region_that_arrived, {&ns[0], &ns[4], 1 << 20}});
    auto blossom2 = ns[3].region_that_arrived->blossom_parent;
    ASSERT_EQ(blossom2->blossom_children.size(), 3);
    blossom2->radius = blossom2->radius + 30;
    mwpm.process_event(RegionHitBoundaryEventData{blossom2, {&ns[0], nullptr, 1 << 30}});
}

TEST(Mwpm, ShatterAndMatchBlossomMatchedToBoundaryMatch) {
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(5, 64)));
    auto& ns = mwpm.flooder.graph.nodes;
    form_nested_blossom_and_match_to_boundary(mwpm);
    auto blossom2 = ns[3].region_that_arrived->blossom_parent;
    auto res = mwpm.shatter_blossom_and_extract_matches(blossom2);
    MatchingResult res_expected(set_bits_to_obs_mask({4, 1, 30}), 1 + 2 + 3 + 4 + 5 + 20 + 30);
    ASSERT_EQ(res, res_expected);
}

TEST(Mwpm, ShatterAndMatchBlossomMatchedToBoundaryMatchEdges) {
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(5, 64)));
    auto& ns = mwpm.flooder.graph.nodes;
    form_nested_blossom_and_match_to_boundary(mwpm);
    auto blossom2 = ns[3].region_that_arrived->blossom_parent;
    std::vector<CompressedEdge> match_edges;
    mwpm.shatter_blossom_and_extract_match_edges(blossom2, match_edges);
    auto obs_expected = set_bits_to_obs_mask({4, 1, 30});
    pm::obs_int obs_mask = 0;
    for (auto& e : match_edges) {
        obs_mask ^= e.obs_mask;
    }
    ASSERT_EQ(obs_mask, obs_expected);
    std::vector<CompressedEdge> expected_match_edges = {
        {&ns[4], &ns[3], 1 << 4}, {&ns[2], &ns[1], 1 << 1}, {&ns[0], nullptr, 1 << 30}};
    ASSERT_EQ(match_edges, expected_match_edges);
}

TEST(Mwpm, MatchingResult) {
    MatchingResult mr;
    ASSERT_EQ(mr.obs_mask, 0);
    ASSERT_EQ(mr.weight, 0);
    MatchingResult mr2(2, 10);
    ASSERT_EQ(mr + mr2, mr2);
    MatchingResult mr3(3, 15);
    ASSERT_EQ(mr2 + mr3, MatchingResult(1, 25));
    MatchingResult mr4(7, 7);
    mr3 += mr4;
    ASSERT_EQ(mr3, MatchingResult(4, 22));
    ASSERT_NE(mr, mr2);
}

TEST(Mwpm, TwoRegionsGrowingThenMatching) {
    Mwpm mwpm(GraphFlooder(MatchingGraph(10, 64)));
    auto& g = mwpm.flooder.graph;
    g.add_boundary_edge(0, 4, {0, 1});
    g.add_edge(0, 1, 100, {0, 2});
    g.add_edge(1, 2, 22, {0, 2});
    g.add_edge(2, 3, 30, {0});
    g.add_edge(3, 4, 10, {0, 3});
    g.add_boundary_edge(4, 1000, {1});
    mwpm.create_detection_event(&mwpm.flooder.graph.nodes[1]);
    mwpm.create_detection_event(&mwpm.flooder.graph.nodes[3]);
    auto e1 = mwpm.flooder.run_until_next_mwpm_notification();
    MwpmEvent e1_expected = RegionHitRegionEventData{
        mwpm.flooder.graph.nodes[1].region_that_arrived,
        mwpm.flooder.graph.nodes[3].region_that_arrived,
        CompressedEdge{&mwpm.flooder.graph.nodes[1], &mwpm.flooder.graph.nodes[3], 4}};
    ASSERT_EQ(e1, e1_expected);
    mwpm.process_event(e1);
    auto e2 = mwpm.flooder.run_until_next_mwpm_notification();
    ASSERT_EQ(e2.event_type, NO_EVENT);
    ASSERT_TRUE(mwpm.flooder.queue.empty());
    std::vector<DetectorNode*> expected_area_1 = {&mwpm.flooder.graph.nodes[1], &mwpm.flooder.graph.nodes[2]};
    ASSERT_EQ(mwpm.flooder.graph.nodes[1].region_that_arrived->shell_area, expected_area_1);
    std::vector<DetectorNode*> expected_area_3 = {&mwpm.flooder.graph.nodes[3], &mwpm.flooder.graph.nodes[4]};
    ASSERT_EQ(mwpm.flooder.graph.nodes[3].region_that_arrived->shell_area, expected_area_3);
    ASSERT_EQ(mwpm.flooder.graph.nodes[1].region_that_arrived->radius, VaryingCT(26 << 2));
    ASSERT_EQ(mwpm.flooder.graph.nodes[3].region_that_arrived->radius, VaryingCT(26 << 2));
}

TEST(Mwpm, RegionHittingMatchThenMatchedToOtherRegion) {
    Mwpm mwpm(GraphFlooder(MatchingGraph(10, 64)));
    auto& g = mwpm.flooder.graph;
    g.add_boundary_edge(0, 1000, {0, 1});
    g.add_edge(0, 1, 8, {0, 2});
    g.add_edge(1, 2, 10, {0, 2});
    g.add_edge(2, 3, 2, {0});
    g.add_edge(3, 4, 4, {0, 3});
    g.add_edge(4, 5, 20, {1});
    g.add_edge(5, 6, 36, {0, 1});
    g.add_boundary_edge(6, 1000, {1});
    mwpm.create_detection_event(&mwpm.flooder.graph.nodes[1]);
    auto r1 = mwpm.flooder.graph.nodes[1].region_that_arrived;
    mwpm.create_detection_event(&mwpm.flooder.graph.nodes[4]);
    auto r4 = mwpm.flooder.graph.nodes[4].region_that_arrived;
    mwpm.create_detection_event(&mwpm.flooder.graph.nodes[5]);
    auto r5 = mwpm.flooder.graph.nodes[5].region_that_arrived;
    mwpm.create_detection_event(&mwpm.flooder.graph.nodes[6]);
    auto r6 = mwpm.flooder.graph.nodes[6].region_that_arrived;
    auto e1 = mwpm.flooder.run_until_next_mwpm_notification();
    MwpmEvent e1_expected = RegionHitRegionEventData{
        r4, r1, CompressedEdge{&mwpm.flooder.graph.nodes[4], &mwpm.flooder.graph.nodes[1], 13}};
    ASSERT_EQ(e1, e1_expected);
    mwpm.process_event(e1);
    auto e2 = mwpm.flooder.run_until_next_mwpm_notification();
    MwpmEvent e2_expected =
        RegionHitRegionEventData{r4, r5, CompressedEdge{&mwpm.flooder.graph.nodes[4], &mwpm.flooder.graph.nodes[5], 2}};
    ASSERT_EQ(e2, e2_expected);
    ASSERT_EQ(mwpm.flooder.queue.cur_time, 12);
    mwpm.process_event(e2);
    auto e3 = mwpm.flooder.run_until_next_mwpm_notification();
    MwpmEvent e3_expected =
        RegionHitRegionEventData{r6, r5, CompressedEdge{&mwpm.flooder.graph.nodes[6], &mwpm.flooder.graph.nodes[5], 3}};
    ASSERT_EQ(e3, e3_expected);
    ASSERT_EQ(mwpm.flooder.queue.cur_time, 18);
    mwpm.process_event(e3);
    auto e4 = mwpm.flooder.run_until_next_mwpm_notification();
    ASSERT_EQ(e4.event_type, NO_EVENT);
    ASSERT_EQ(r1->radius, VaryingCT(14 << 2));
    std::vector<DetectorNode*> area_1 = {
        &mwpm.flooder.graph.nodes[1],
        &mwpm.flooder.graph.nodes[0],
        &mwpm.flooder.graph.nodes[2],
        &mwpm.flooder.graph.nodes[3]};
    ASSERT_EQ(r1->shell_area, area_1);
    ASSERT_EQ(r4->radius, VaryingCT(2 << 2));
    std::vector<DetectorNode*> area_4 = {&mwpm.flooder.graph.nodes[4]};
    ASSERT_EQ(r4->shell_area, area_4);
    for (auto ri : {r1, r4, r5, r6})
        ASSERT_EQ(ri->alt_tree_node, nullptr);
    ASSERT_EQ(r1->match, (Match{r4, CompressedEdge{&mwpm.flooder.graph.nodes[1], &mwpm.flooder.graph.nodes[4], 13}}));
    ASSERT_EQ(r4->match, (Match{r1, CompressedEdge{&mwpm.flooder.graph.nodes[4], &mwpm.flooder.graph.nodes[1], 13}}));
    ASSERT_EQ(r5->match, (Match{r6, CompressedEdge{&mwpm.flooder.graph.nodes[5], &mwpm.flooder.graph.nodes[6], 3}}));
    ASSERT_EQ(r6->match, (Match{r5, CompressedEdge{&mwpm.flooder.graph.nodes[6], &mwpm.flooder.graph.nodes[5], 3}}));
}

TEST(Mwpm, RegionHittingMatchFormingBlossomThenMatchingToBoundary) {
    size_t num_nodes = 100;
    Mwpm mwpm{GraphFlooder(MatchingGraph(num_nodes, 64))};
    auto& g = mwpm.flooder.graph;
    g.add_boundary_edge(0, 2, {0});
    for (size_t i = 0; i < num_nodes - 1; i++)
        g.add_edge(i, i + 1, 2, obs_mask_to_set_bits(i));
    g.add_boundary_edge(num_nodes - 1, 2, {1});
    mwpm.create_detection_event(&mwpm.flooder.graph.nodes[40]);
    auto r40 = mwpm.flooder.graph.nodes[40].region_that_arrived;
    mwpm.create_detection_event(&mwpm.flooder.graph.nodes[42]);
    auto r42 = mwpm.flooder.graph.nodes[42].region_that_arrived;
    mwpm.create_detection_event(&mwpm.flooder.graph.nodes[50]);
    auto r50 = mwpm.flooder.graph.nodes[50].region_that_arrived;
    mwpm.process_event(mwpm.flooder.run_until_next_mwpm_notification());
    mwpm.process_event(mwpm.flooder.run_until_next_mwpm_notification());
    mwpm.process_event(mwpm.flooder.run_until_next_mwpm_notification());
    mwpm.process_event(mwpm.flooder.run_until_next_mwpm_notification());
    ASSERT_EQ(mwpm.flooder.queue.cur_time, 94);
    std::vector<RegionEdge> expected_blossom_regions = {
        {r40, {&mwpm.flooder.graph.nodes[40], &mwpm.flooder.graph.nodes[42], 40 ^ 41}},
        {r42, {&mwpm.flooder.graph.nodes[42], &mwpm.flooder.graph.nodes[50], 42 ^ 43 ^ 44 ^ 45 ^ 46 ^ 47 ^ 48 ^ 49}},
        {r50,
         {&mwpm.flooder.graph.nodes[50],
          &mwpm.flooder.graph.nodes[40],
          40 ^ 41 ^ 42 ^ 43 ^ 44 ^ 45 ^ 46 ^ 47 ^ 48 ^ 49}},

    };
    auto blossom = mwpm.flooder.graph.nodes[40].region_that_arrived->blossom_parent;
    auto actual_blossom_regions = blossom->blossom_children;
    ASSERT_EQ(actual_blossom_regions, expected_blossom_regions);
    obs_int obs = 0;
    for (int i = 0; i < 40; i++)
        obs ^= i;
    ASSERT_EQ(blossom->match, (Match{nullptr, {&mwpm.flooder.graph.nodes[40], nullptr, obs ^ 1}}));
}

TEST(GraphFlooder, CreateRegion) {
    Mwpm mwpm{GraphFlooder(MatchingGraph(5, 64))};
    auto& flooder = mwpm.flooder;
    auto& g = flooder.graph;
    g.add_boundary_edge(0, 3, {});
    g.add_edge(0, 1, 5, {});
    g.add_edge(1, 2, 11, {});
    g.add_edge(2, 3, 100, {});
    g.add_boundary_edge(3, 1000, {});
    mwpm.create_detection_event(&flooder.graph.nodes[0]);
    mwpm.create_detection_event(&flooder.graph.nodes[2]);
    mwpm.create_detection_event(&flooder.graph.nodes[3]);
    ASSERT_EQ(flooder.queue.dequeue(), FloodCheckEvent(&flooder.graph.nodes[0], cyclic_time_int{3}));
    ASSERT_EQ(flooder.queue.dequeue(), FloodCheckEvent(&flooder.graph.nodes[2], cyclic_time_int{11}));
    ASSERT_EQ(flooder.graph.nodes[0].region_that_arrived->shell_area[0], &flooder.graph.nodes[0]);
    ASSERT_EQ(flooder.graph.nodes[0].radius_of_arrival, 0);
    ASSERT_EQ(flooder.graph.nodes[0].reached_from_source, &flooder.graph.nodes[0]);
    ASSERT_EQ(flooder.queue.dequeue(), FloodCheckEvent(&flooder.graph.nodes[3], cyclic_time_int{50}));
    ASSERT_TRUE(flooder.queue.empty());
}

TEST(GraphFlooder, RegionGrowingToBoundary) {
    Mwpm mwpm{GraphFlooder(MatchingGraph(10, 64))};
    auto& flooder = mwpm.flooder;
    auto& g = flooder.graph;
    g.add_boundary_edge(0, 2, {0, 1});
    g.add_edge(0, 1, 10, {0, 2});
    g.add_edge(1, 2, 21, {});
    g.add_edge(2, 3, 100, {0});
    g.add_edge(3, 4, 7, {0, 3});
    g.add_boundary_edge(4, 5, {1});
    mwpm.create_detection_event(&flooder.graph.nodes[2]);
    ASSERT_EQ(
        flooder.run_until_next_mwpm_notification(),
        (MwpmEvent{
            RegionHitBoundaryEventData{
                flooder.graph.nodes[2].region_that_arrived,
                CompressedEdge{&flooder.graph.nodes[2], nullptr, 6},
            },
        }));

    flooder.queue.dequeue();
    auto t1 = flooder.queue.dequeue();
    ASSERT_EQ(t1, FloodCheckEvent(&flooder.graph.nodes[2], cyclic_time_int{100}));
    flooder.queue.enqueue(t1);
    auto e2 = flooder.run_until_next_mwpm_notification();
    MwpmEvent e2_exp = RegionHitBoundaryEventData{
        flooder.graph.nodes[2].region_that_arrived,
        CompressedEdge{&flooder.graph.nodes[2], nullptr, 10},
    };
    ASSERT_EQ(e2, e2_exp);
    flooder.queue.dequeue();
    ASSERT_EQ(flooder.queue.dequeue(), FloodCheckEvent(cyclic_time_int{0}));
}

TEST(GraphFlooder, RegionHitRegion) {
    Mwpm mwpm{GraphFlooder(MatchingGraph(10, 64))};
    auto& flooder = mwpm.flooder;
    auto& g = flooder.graph;
    g.add_boundary_edge(0, 2000, {0, 1});
    g.add_edge(0, 1, 10, {0, 2});
    g.add_edge(1, 2, 24, {});
    g.add_edge(2, 3, 38, {0});
    g.add_edge(3, 4, 26, {0, 3});
    g.add_boundary_edge(4, 1000, {1});
    mwpm.create_detection_event(&flooder.graph.nodes[2]);
    mwpm.create_detection_event(&flooder.graph.nodes[4]);
    auto e1 = flooder.run_until_next_mwpm_notification();
    ASSERT_EQ(flooder.queue.cur_time, 32);
    MwpmEvent e1_exp = RegionHitRegionEventData{
        flooder.graph.nodes[4].region_that_arrived,
        flooder.graph.nodes[2].region_that_arrived,
        CompressedEdge{&flooder.graph.nodes[4], &flooder.graph.nodes[2], 8}};
    ASSERT_EQ(e1, e1_exp);
}

TEST(GraphFlooder, RegionGrowingThenFrozenThenStartShrinking) {
    Mwpm mwpm{GraphFlooder(MatchingGraph(10, 64))};
    auto& flooder = mwpm.flooder;
    auto& g = flooder.graph;
    g.add_boundary_edge(0, 4, {0, 1});
    g.add_edge(0, 1, 10, {0, 2});
    g.add_edge(1, 2, 22, {});
    g.add_edge(2, 3, 30, {0});
    g.add_edge(3, 4, 50, {0, 3});
    g.add_boundary_edge(4, 100, {1});
    mwpm.create_detection_event(&flooder.graph.nodes[2]);
    auto e1 = flooder.run_until_next_mwpm_notification();
    ASSERT_EQ(
        e1,
        MwpmEvent(RegionHitBoundaryEventData{
            flooder.graph.nodes[2].region_that_arrived,
            CompressedEdge{&flooder.graph.nodes[2], nullptr, 6},
        }));
    ASSERT_EQ(flooder.queue.cur_time, 36);
    flooder.set_region_frozen(*flooder.graph.nodes[2].region_that_arrived);
    auto e2 = flooder.run_until_next_mwpm_notification();
    ASSERT_EQ(e2.event_type, NO_EVENT);
    ASSERT_EQ(flooder.queue.cur_time, 80);
    ASSERT_TRUE(flooder.queue.empty());

    std::vector<DetectorNode*> expected_area = {
        &flooder.graph.nodes[2], &flooder.graph.nodes[1], &flooder.graph.nodes[3], &flooder.graph.nodes[0]};
    ASSERT_EQ(flooder.graph.nodes[2].region_that_arrived->shell_area, expected_area);
    flooder.set_region_shrinking(*flooder.graph.nodes[2].region_that_arrived);
    ASSERT_EQ(
        flooder.dequeue_valid(), FloodCheckEvent(flooder.graph.nodes[2].region_that_arrived, cyclic_time_int{80 + 4}));
}
