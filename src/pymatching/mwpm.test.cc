#include "pymatching/alternating_tree.h"
#include "pymatching/graph_fill_region.h"
#include "pymatching/mwpm.h"

#include <gtest/gtest.h>

using namespace pm;

struct MwpmAltTreeTestData {
    std::vector<DetectorNode>* nodes;
    explicit MwpmAltTreeTestData(std::vector<DetectorNode>* nodes);
    AltTreeEdge t(
        std::vector<AltTreeEdge> children, size_t inner_region_id, size_t outer_region_id, bool root = false) const;
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
    return RegionHitRegionEventData{
        ns[i].region_that_arrived, ns[j].region_that_arrived, {&ns[i], &ns[j], obs_mask}};
}

MwpmAltTreeTestData::MwpmAltTreeTestData(std::vector<DetectorNode>* nodes) : nodes(nodes) {
}

AltTreeEdge MwpmAltTreeTestData::t(
    std::vector<AltTreeEdge> children, size_t inner_region_id, size_t outer_region_id, bool root) const {
    AltTreeNode* node;
    CompressedEdge parent_ce(nullptr, nullptr, 0);
    if (root) {
        node = new AltTreeNode((*nodes)[outer_region_id].region_that_arrived);
    } else {
        node = new AltTreeNode(
            (*nodes)[inner_region_id].region_that_arrived,
            (*nodes)[outer_region_id].region_that_arrived,
            {&(*nodes)[inner_region_id], &(*nodes)[outer_region_id], 0});
    }
    auto edge = AltTreeEdge(node, parent_ce);
    for (auto& child : children) {
        child.edge.loc_from = &(*nodes)[outer_region_id];
        child.edge.loc_to = child.alt_tree_node->inner_to_outer_edge.loc_from;
        edge.alt_tree_node->add_child(child);
    }
    return edge;
}

TEST(Mwpm, BlossomCreatedThenShattered) {
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(10)));
    auto &g = mwpm.flooder.graph;
    g.add_edge(0, 1, 10, 1);
    g.add_edge(1, 4, 20, 2);
    g.add_edge(4, 3, 20, 3);
    g.add_edge(3, 2, 12, 4);
    g.add_edge(0, 2, 16, 5);
    g.add_edge(4, 5, 50, 6);
    g.add_edge(2, 6, 100, 7);
    g.add_boundary_edge(5, 36, 8);
    for (size_t i = 0; i < 7; i++) {
        mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[i]);
    }
    auto& ns = mwpm.flooder.graph.nodes;

    auto dequeued_and_processed = [&](size_t expected_time){
        auto ev = mwpm.flooder.run_until_next_mwpm_notification();
        EXPECT_EQ(expected_time, mwpm.flooder.queue.cur_time);
        mwpm.process_event(ev);
        return ev;
    };

    // 0:1 and 2:3 pair up.
    ASSERT_EQ(
        dequeued_and_processed(5),
        MwpmEvent(RegionHitRegionEventData{
            ns[0].region_that_arrived,
            ns[1].region_that_arrived,
            {&ns[0], &ns[1], 1}}));
    ASSERT_EQ(
        dequeued_and_processed(6),
        MwpmEvent(RegionHitRegionEventData{
            ns[2].region_that_arrived,
            ns[3].region_that_arrived,
            {&ns[2], &ns[3], 4}}));

    // 4 hits frozen 3 creating 4-3-2 tree.
    ASSERT_EQ(
        dequeued_and_processed(14),
        MwpmEvent(RegionHitRegionEventData{
            ns[3].region_that_arrived,
            ns[4].region_that_arrived,
            {&ns[3], &ns[4], 3}}));
    // 4 hits frozen 1 creating 0-1-4-3-2 tree.
    ASSERT_EQ(
        dequeued_and_processed(15),
        MwpmEvent(RegionHitRegionEventData{
            ns[1].region_that_arrived,
            ns[4].region_that_arrived,
            {&ns[1], &ns[4], 2}}));

    // 0 hits 2 creating blossom (0-1-4-3-2).
    ASSERT_EQ(
        dequeued_and_processed(17),
        MwpmEvent(RegionHitRegionEventData{
            ns[0].region_that_arrived,
            ns[2].region_that_arrived,
            {&ns[0], &ns[2], 5}}));
    auto *blossom = ns[0].region_that_arrived->blossom_parent;
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
        MwpmEvent(RegionHitRegionEventData{
            blossom,
            ns[5].region_that_arrived,
            {&ns[4], &ns[5], 6}}));
    ASSERT_EQ(blossom->radius, Varying(8 << 2));

    // Region 6 collides with matched blossom
    ASSERT_EQ(
        dequeued_and_processed(83),
        MwpmEvent(RegionHitRegionEventData{
            blossom,
            ns[6].region_that_arrived,
            {&ns[2], &ns[6], 7}}));
    ASSERT_EQ(
        ns[6].region_that_arrived->alt_tree_node->children[0].alt_tree_node->inner_region,
        blossom);

    // Blossom shatters
    ASSERT_EQ(
        dequeued_and_processed(91),
        MwpmEvent(BlossomShatterEventData{
            blossom,
            ns[2].region_that_arrived,
            ns[4].region_that_arrived}));
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
    ASSERT_EQ(ns[5].region_that_arrived->match, Match(nullptr, {&ns[5], nullptr, 8}));
}

TEST(Mwpm, BlossomShatterDrivenWithoutFlooder) {
    size_t n = 10;
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(n + 3)));
    auto& ns = mwpm.flooder.graph.nodes;
    for (size_t i = 0; i < n + 3; i++)
        mwpm.flooder.create_region(&ns[i]);

    // Pair up
    for (size_t i = 0; i < n; i += 2)
        mwpm.process_event(rhr(ns, i, i + 1));
    // Form an alternating path
    for (size_t i = 1; i < n; i += 2)
        mwpm.process_event(rhr(ns, n - i, n - i + 1));
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

    MwpmAltTreeTestData d(&ns);
    auto expected_tree = d.t({d.t({d.t({d.t({d.t({}, 5, 11)}, 3, 4)}, 1, 2)}, 10, 0)}, -1, 12, true).alt_tree_node;
    ASSERT_EQ(*actual_tree, *expected_tree);
}

TEST(Mwpm, BranchingTreeFormsBlossomThenHitsBoundaryMatch) {
    size_t n = 14;
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(n)));
    auto& ns = mwpm.flooder.graph.nodes;
    for (size_t i = 0; i < n; i++)
        mwpm.flooder.create_region(&ns[i]);

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
    MwpmAltTreeTestData d(&ns);
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
    mwpm.process_event(RegionHitBoundaryEventData{ns[13].region_that_arrived,
                                                  {&ns[13], nullptr, 0}});
    // Blossom matches to region 13
    mwpm.process_event(
        MwpmEvent(RegionHitRegionEventData{blossom, ns[13].region_that_arrived, {&ns[4], &ns[13], 1}}));
    ASSERT_TRUE(regions_matched(ns, 5, 9));
    ASSERT_TRUE(regions_matched(ns, 8, 12));
    ASSERT_EQ(blossom->match, Match(ns[13].region_that_arrived, {&ns[4], &ns[13], 1}));
    ASSERT_EQ(ns[13].region_that_arrived->match, Match(blossom, {&ns[13], &ns[4], 1}));
}

TEST(Mwpm, BoundaryMatchHitsTree) {
    size_t n = 6;
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(n)));
    auto& ns = mwpm.flooder.graph.nodes;
    for (size_t i = 0; i < n; i++)
        mwpm.flooder.create_region(&ns[i]);
    // Form alternating tree
    mwpm.process_event(rhr(ns, 1, 2));
    mwpm.process_event(rhr(ns, 3, 4));
    mwpm.process_event(rhr(ns, 0, 1));
    mwpm.process_event(rhr(ns, 0, 3));
    // 5 matches to boundary
    mwpm.process_event(RegionHitBoundaryEventData{ns[5].region_that_arrived,
                                                  {&ns[5], nullptr, 1}});
    // Tree matches to 5
    mwpm.process_event(rhr(ns, 5, 2));
    ASSERT_TRUE(regions_matched(ns, 5, 2));
    ASSERT_TRUE(regions_matched(ns, 1, 0));
    ASSERT_TRUE(regions_matched(ns, 3, 4));
};

TEST(Mwpm, ShatterBlossomAndExtractMatchesForPair) {
    size_t num_nodes = 20;
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(num_nodes)));
    auto &g = mwpm.flooder.graph;
    for (size_t i = 0; i < num_nodes - 1; i++)
        g.add_edge(i, i + 1, 2, i);
    auto& ns = mwpm.flooder.graph.nodes;
    mwpm.flooder.create_region(&ns[6]);
    mwpm.flooder.create_region(&ns[13]);
    auto e = mwpm.flooder.run_until_next_mwpm_notification();
    mwpm.process_event(e);
    ASSERT_EQ(ns[6].region_that_arrived->match.region, ns[13].region_that_arrived);
    auto res = mwpm.shatter_blossom_and_extract_matches(ns[6].region_that_arrived);
    ASSERT_EQ(res, MatchingResult(6 ^ 7 ^ 8 ^ 9 ^ 10 ^ 11 ^ 12, 14));
    for (auto& n : ns) {
        ASSERT_EQ(n.region_that_arrived, nullptr);
        ASSERT_EQ(n.distance_from_source, 0);
        ASSERT_EQ(n.observables_crossed_from_source, 0);
        ASSERT_EQ(n.reached_from_source, nullptr);
    }
}

obs_int set_bits_to_obs_mask(const std::vector<int>& set_bits) {
    obs_int obs_mask = 0;
    for (auto i : set_bits)
        obs_mask ^= (1 << i);
    return obs_mask;
};

TEST(Mwpm, ShatterBlossomAndExtractMatches) {
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(12)));

    auto& ns = mwpm.flooder.graph.nodes;
    cumulative_time_int w = 0;
    for (auto& n : ns) {
        w += 1;
        mwpm.flooder.create_region(&n);
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
    auto res = mwpm.shatter_blossom_and_extract_matches(blossom1);

    ASSERT_EQ(res.obs_mask, set_bits_to_obs_mask({3, 4, 15, 6, 14, 9}));
    ASSERT_EQ(res.weight, 2 + 3 + 4 + 5 + 1 + 2 + 6 + 3 + 11 + 7 + 8 + 9 + 10 + 5 + 11 + 12);
}

TEST(Mwpm, ShatterAndMatchBlossomMatchedToBoundary) {
    auto mwpm = Mwpm(GraphFlooder(MatchingGraph(5)));
    auto& ns = mwpm.flooder.graph.nodes;
    cumulative_time_int w = 0;
    for (auto& n : ns) {
        w += 1;
        mwpm.flooder.create_region(&n);
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
    mwpm.process_event(RegionHitBoundaryEventData{blossom2,
                                                  {&ns[0], nullptr, 1 << 30}});
    auto res = mwpm.shatter_blossom_and_extract_matches(blossom2);
    MatchingResult res_expected(set_bits_to_obs_mask({4, 1, 30}), 1 + 2 + 3 + 4 + 5 + 20 + 30);
    ASSERT_EQ(res, res_expected);
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
