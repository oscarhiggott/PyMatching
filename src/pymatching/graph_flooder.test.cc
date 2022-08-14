#include "pymatching/graph_flooder.h"

#include "gtest/gtest.h"

#include "pymatching/events.h"
#include "pymatching/mwpm.h"

TEST(GraphFlooder, PriorityQueue) {
    pm::GraphFlooder flooder(pm::MatchingGraph(10), 1024);
    auto &graph = flooder.graph;
    graph.add_edge(0, 1, 10, 0);
    graph.add_edge(1, 2, 10, 0);
    graph.add_edge(2, 3, 10, 0);
    graph.add_edge(3, 4, 10, 0);
    graph.add_edge(4, 5, 10, 0);
    graph.add_edge(5, 0, 10, 0);

    flooder.queue.enqueue(pm::TentativeEvent(pm::TentativeNeighborInteractionEventData{
        &graph.nodes[0], 0, &graph.nodes[1], 1
    }, 10, 0x0));
    flooder.queue.enqueue(pm::TentativeEvent(pm::TentativeNeighborInteractionEventData{
        &graph.nodes[1], 0, &graph.nodes[2], 1
    }, 8, 0x0));
    flooder.queue.enqueue(pm::TentativeEvent(pm::TentativeNeighborInteractionEventData{
        &graph.nodes[2], 0, &graph.nodes[3], 1
    }, 5, 0x0));

    pm::GraphFillRegion gfr;
    flooder.queue.enqueue(pm::TentativeEvent(pm::TentativeRegionShrinkEventData{&gfr}, 70, 0x0));

    flooder.queue.enqueue(pm::TentativeEvent(pm::TentativeNeighborInteractionEventData{
        &graph.nodes[4], 0, &graph.nodes[5], 1
    }, 100, 0x0));
    auto e = flooder.queue.dequeue_valid();
    ASSERT_EQ(e.time, 5);
    ASSERT_EQ(e.tentative_event_type, pm::INTERACTION);
    ASSERT_EQ(e.neighbor_interaction_event_data.detector_node_1, &graph.nodes[2]);
    ASSERT_EQ(flooder.queue.dequeue_valid().time, 8);
    ASSERT_EQ(flooder.queue.dequeue_valid().time, 10);
    auto e6 = flooder.queue.dequeue_valid();
    ASSERT_EQ(e6.time, 70);
    ASSERT_EQ(e6.tentative_event_type, pm::SHRINKING);
    ASSERT_EQ(e6.region_shrink_event_data.region, &gfr);
    ASSERT_TRUE(e6.is_still_valid());
    ASSERT_EQ(flooder.queue.dequeue_valid().time, 100);
}

TEST(GraphFlooder, CreateRegion) {
    pm::GraphFlooder flooder(pm::MatchingGraph(5), 2048);
    auto &g = flooder.graph;
    g.add_boundary_edge(0, 3, 0);
    g.add_edge(0, 1, 5, 0);
    g.add_edge(1, 2, 11, 0);
    g.add_edge(2, 3, 100, 0);
    g.add_boundary_edge(3, 1000, 0);
    flooder.create_region(&flooder.graph.nodes[0]);
    flooder.create_region(&flooder.graph.nodes[2]);
    flooder.create_region(&flooder.graph.nodes[3]);
    ASSERT_EQ(flooder.queue.dequeue_valid(), pm::TentativeEvent(pm::TentativeNeighborInteractionEventData{
        &flooder.graph.nodes[0], 0, nullptr, (size_t)-1
    }, 3, 0x0));
    ASSERT_EQ(flooder.queue.dequeue_valid(), pm::TentativeEvent(pm::TentativeNeighborInteractionEventData{
        &flooder.graph.nodes[0], 1, &flooder.graph.nodes[1], 0
    }, 5, 0x1));
    ASSERT_EQ(flooder.graph.nodes[0].region_that_arrived->shell_area[0], &flooder.graph.nodes[0]);
    ASSERT_EQ(flooder.graph.nodes[0].distance_from_source, 0);
    ASSERT_EQ(flooder.graph.nodes[0].reached_from_source, &flooder.graph.nodes[0]);
    ASSERT_EQ(flooder.queue.dequeue_valid(), pm::TentativeEvent(pm::TentativeNeighborInteractionEventData{
        &flooder.graph.nodes[2], 0, &flooder.graph.nodes[1], 1
    }, 11, 0x2));
    ASSERT_EQ(flooder.queue.dequeue_valid(), pm::TentativeEvent(pm::TentativeNeighborInteractionEventData{
        &flooder.graph.nodes[3], 0, &flooder.graph.nodes[2], 1
    }, 50, 0x4));
    // t=100 event skipped over because it was invalidated.
    ASSERT_EQ(flooder.queue.dequeue_valid(), pm::TentativeEvent(pm::TentativeNeighborInteractionEventData{
        &flooder.graph.nodes[3], 1, nullptr, (size_t)-1
    }, 1000, 0x5));
    ASSERT_TRUE(flooder.queue.empty());
}

TEST(GraphFlooder, RegionGrowingToBoundary) {
    pm::GraphFlooder flooder(pm::MatchingGraph(10), 1024);
    auto &g = flooder.graph;
    g.add_boundary_edge(0, 2, 3);
    g.add_edge(0, 1, 10, 5);
    g.add_edge(1, 2, 21, 0);
    g.add_edge(2, 3, 100, 1);
    g.add_edge(3, 4, 7, 9);
    g.add_boundary_edge(4, 5, 2);
    flooder.create_region(&flooder.graph.nodes[2]);
    auto e1 = flooder.next_event();
    pm::MwpmEvent e1_exp = pm::RegionHitBoundaryEventData(
        flooder.graph.nodes[2].region_that_arrived, pm::CompressedEdge(&flooder.graph.nodes[2], nullptr, 6));
    ASSERT_EQ(e1, e1_exp);
    auto t1 = flooder.queue.dequeue_valid();
    ASSERT_EQ(t1, pm::TentativeEvent(pm::TentativeNeighborInteractionEventData{
        &flooder.graph.nodes[2], 1, &flooder.graph.nodes[3], 0
    }, 100, 0x1));
    flooder.queue.enqueue(pm::TentativeEvent(t1));
    auto e2 = flooder.next_event();
    pm::MwpmEvent e2_exp = pm::RegionHitBoundaryEventData(
        flooder.graph.nodes[2].region_that_arrived, pm::CompressedEdge(&flooder.graph.nodes[2], nullptr, 10));
    ASSERT_EQ(e2, e2_exp);
    auto e3 = flooder.next_event();
    ASSERT_EQ(e3.event_type, pm::NO_EVENT);
}

TEST(GraphFlooder, RegionHitRegion) {
    pm::GraphFlooder flooder(pm::MatchingGraph(10), 2048);
    auto &g = flooder.graph;
    g.add_boundary_edge(0, 2000, 3);
    g.add_edge(0, 1, 10, 5);
    g.add_edge(1, 2, 24, 0);
    g.add_edge(2, 3, 38, 1);
    g.add_edge(3, 4, 26, 9);
    g.add_boundary_edge(4, 1000, 2);
    flooder.create_region(&flooder.graph.nodes[2]);
    flooder.create_region(&flooder.graph.nodes[4]);
    auto e1 = flooder.next_event();
    ASSERT_EQ(flooder.queue.cur_time, 32);
    pm::MwpmEvent e1_exp = pm::RegionHitRegionEventData(
        flooder.graph.nodes[4].region_that_arrived,
        flooder.graph.nodes[2].region_that_arrived,
        pm::CompressedEdge(&flooder.graph.nodes[4], &flooder.graph.nodes[2], 8));
    ASSERT_EQ(e1, e1_exp);
}

TEST(GraphFlooder, RegionGrowingThenFrozenThenStartShrinking) {
    pm::GraphFlooder flooder(pm::MatchingGraph(10), 1024);
    auto &g = flooder.graph;
    g.add_boundary_edge(0, 4, 3);
    g.add_edge(0, 1, 10, 5);
    g.add_edge(1, 2, 22, 0);
    g.add_edge(2, 3, 30, 1);
    g.add_edge(3, 4, 50, 9);
    g.add_boundary_edge(4, 100, 2);
    flooder.create_region(&flooder.graph.nodes[2]);
    auto e1 = flooder.next_event();
    ASSERT_EQ(
        e1,
        pm::MwpmEvent(pm::RegionHitBoundaryEventData(
            flooder.graph.nodes[2].region_that_arrived, pm::CompressedEdge(&flooder.graph.nodes[2], nullptr, 6))));
    ASSERT_EQ(flooder.queue.cur_time, 36);
    flooder.set_region_frozen(*flooder.graph.nodes[2].region_that_arrived);
    auto e2 = flooder.next_event();
    ASSERT_EQ(e2.event_type, pm::NO_EVENT);
    ASSERT_EQ(flooder.queue.cur_time, 80);
    ASSERT_TRUE(flooder.queue.empty());

    std::vector<pm::DetectorNode*> expected_area = {
        &flooder.graph.nodes[2], &flooder.graph.nodes[1], &flooder.graph.nodes[3], &flooder.graph.nodes[0]};
    ASSERT_EQ(flooder.graph.nodes[2].region_that_arrived->shell_area, expected_area);
    flooder.set_region_shrinking(*flooder.graph.nodes[2].region_that_arrived);
    ASSERT_EQ(flooder.queue.dequeue_valid(), pm::TentativeEvent(pm::TentativeRegionShrinkEventData{flooder.graph.nodes[2].region_that_arrived}, 80 + 4, 0x5));
}

TEST(GraphFlooder, TwoRegionsGrowingThenMatching) {
    pm::Mwpm mwpm(pm::GraphFlooder(pm::MatchingGraph(10), 2048));
    auto &g = mwpm.flooder.graph;
    g.add_boundary_edge(0, 4, 3);
    g.add_edge(0, 1, 100, 5);
    g.add_edge(1, 2, 22, 5);
    g.add_edge(2, 3, 30, 1);
    g.add_edge(3, 4, 10, 9);
    g.add_boundary_edge(4, 1000, 2);
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[1]);
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[3]);
    auto e1 = mwpm.flooder.next_event();
    pm::MwpmEvent e1_expected = pm::RegionHitRegionEventData(
        mwpm.flooder.graph.nodes[1].region_that_arrived,
        mwpm.flooder.graph.nodes[3].region_that_arrived,
        pm::CompressedEdge(&mwpm.flooder.graph.nodes[1], &mwpm.flooder.graph.nodes[3], 4));
    ASSERT_EQ(e1, e1_expected);
    mwpm.process_event(e1);
    auto e2 = mwpm.flooder.next_event();
    ASSERT_EQ(e2.event_type, pm::NO_EVENT);
    ASSERT_TRUE(mwpm.flooder.queue.empty());
    std::vector<pm::DetectorNode*> expected_area_1 = {&mwpm.flooder.graph.nodes[1], &mwpm.flooder.graph.nodes[2]};
    ASSERT_EQ(mwpm.flooder.graph.nodes[1].region_that_arrived->shell_area, expected_area_1);
    std::vector<pm::DetectorNode*> expected_area_3 = {&mwpm.flooder.graph.nodes[3], &mwpm.flooder.graph.nodes[4]};
    ASSERT_EQ(mwpm.flooder.graph.nodes[3].region_that_arrived->shell_area, expected_area_3);
    ASSERT_EQ(mwpm.flooder.graph.nodes[1].region_that_arrived->radius, pm::Varying(26 << 2));
    ASSERT_EQ(mwpm.flooder.graph.nodes[3].region_that_arrived->radius, pm::Varying(26 << 2));
}

TEST(GraphFlooder, RegionHittingMatchThenMatchedToOtherRegion) {
    pm::Mwpm mwpm(pm::GraphFlooder(pm::MatchingGraph(10), 2048));
    auto &g = mwpm.flooder.graph;
    g.add_boundary_edge(0, 1000, 3);
    g.add_edge(0, 1, 8, 5);
    g.add_edge(1, 2, 10, 5);
    g.add_edge(2, 3, 2, 1);
    g.add_edge(3, 4, 4, 9);
    g.add_edge(4, 5, 20, 2);
    g.add_edge(5, 6, 36, 3);
    g.add_boundary_edge(6, 1000, 2);
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[1]);
    auto r1 = mwpm.flooder.graph.nodes[1].region_that_arrived;
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[4]);
    auto r4 = mwpm.flooder.graph.nodes[4].region_that_arrived;
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[5]);
    auto r5 = mwpm.flooder.graph.nodes[5].region_that_arrived;
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[6]);
    auto r6 = mwpm.flooder.graph.nodes[6].region_that_arrived;
    auto e1 = mwpm.flooder.next_event();
    pm::MwpmEvent e1_expected = pm::RegionHitRegionEventData(
        r4, r1, pm::CompressedEdge(&mwpm.flooder.graph.nodes[4], &mwpm.flooder.graph.nodes[1], 13));
    ASSERT_EQ(e1, e1_expected);
    mwpm.process_event(e1);
    auto e2 = mwpm.flooder.next_event();
    pm::MwpmEvent e2_expected = pm::RegionHitRegionEventData(
        r4, r5, pm::CompressedEdge(&mwpm.flooder.graph.nodes[4], &mwpm.flooder.graph.nodes[5], 2));
    ASSERT_EQ(e2, e2_expected);
    ASSERT_EQ(mwpm.flooder.queue.cur_time, 12);
    mwpm.process_event(e2);
    auto e3 = mwpm.flooder.next_event();
    pm::MwpmEvent e3_expected = pm::RegionHitRegionEventData(
        r6, r5, pm::CompressedEdge(&mwpm.flooder.graph.nodes[6], &mwpm.flooder.graph.nodes[5], 3));
    ASSERT_EQ(e3, e3_expected);
    ASSERT_EQ(mwpm.flooder.queue.cur_time, 18);
    mwpm.process_event(e3);
    auto e4 = mwpm.flooder.next_event();
    ASSERT_EQ(e4.event_type, pm::NO_EVENT);
    ASSERT_EQ(r1->radius, pm::Varying(14 << 2));
    std::vector<pm::DetectorNode*> area_1 = {
        &mwpm.flooder.graph.nodes[1],
        &mwpm.flooder.graph.nodes[0],
        &mwpm.flooder.graph.nodes[2],
        &mwpm.flooder.graph.nodes[3]};
    ASSERT_EQ(r1->shell_area, area_1);
    ASSERT_EQ(r4->radius, pm::Varying(2 << 2));
    std::vector<pm::DetectorNode*> area_4 = {&mwpm.flooder.graph.nodes[4]};
    ASSERT_EQ(r4->shell_area, area_4);
    for (auto ri : {r1, r4, r5, r6})
        ASSERT_EQ(ri->alt_tree_node, nullptr);
    ASSERT_EQ(
        r1->match, pm::Match(r4, pm::CompressedEdge(&mwpm.flooder.graph.nodes[1], &mwpm.flooder.graph.nodes[4], 13)));
    ASSERT_EQ(
        r4->match, pm::Match(r1, pm::CompressedEdge(&mwpm.flooder.graph.nodes[4], &mwpm.flooder.graph.nodes[1], 13)));
    ASSERT_EQ(
        r5->match, pm::Match(r6, pm::CompressedEdge(&mwpm.flooder.graph.nodes[5], &mwpm.flooder.graph.nodes[6], 3)));
    ASSERT_EQ(
        r6->match, pm::Match(r5, pm::CompressedEdge(&mwpm.flooder.graph.nodes[6], &mwpm.flooder.graph.nodes[5], 3)));
}

TEST(GraphFlooder, RegionHittingMatchFormingBlossomThenMatchingToBoundary) {
    size_t num_nodes = 100;
    pm::Mwpm mwpm{pm::GraphFlooder(pm::MatchingGraph(num_nodes), 2048)};
    auto &g = mwpm.flooder.graph;
    g.add_boundary_edge(0, 2, 1);
    for (size_t i = 0; i < num_nodes - 1; i++)
        g.add_edge(i, i + 1, 2, i);
    g.add_boundary_edge(num_nodes - 1, 2, 2);
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[40]);
    auto r40 = mwpm.flooder.graph.nodes[40].region_that_arrived;
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[42]);
    auto r42 = mwpm.flooder.graph.nodes[42].region_that_arrived;
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[50]);
    auto r50 = mwpm.flooder.graph.nodes[50].region_that_arrived;
    mwpm.process_event(mwpm.flooder.next_event());
    mwpm.process_event(mwpm.flooder.next_event());
    mwpm.process_event(mwpm.flooder.next_event());
    mwpm.process_event(mwpm.flooder.next_event());
    ASSERT_EQ(mwpm.flooder.queue.cur_time, 94);
    std::vector<pm::RegionEdge> expected_blossom_regions = {
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
    pm::obs_int obs = 0;
    for (int i = 0; i < 40; i++)
        obs ^= i;
    ASSERT_EQ(blossom->match, pm::Match(nullptr, {&mwpm.flooder.graph.nodes[40], nullptr, obs ^ 1}));
}
