#include "pymatching/fill_match/flooder/graph_flooder.h"

#include "gtest/gtest.h"

#include "pymatching/fill_match/flooder/graph_fill_region.h"
#include "pymatching/fill_match/ints.h"
#include "pymatching/fill_match/flooder_matcher_interop/mwpm_event.h"

using namespace pm;

TEST(GraphFlooder, PriorityQueue) {
    GraphFlooder flooder(MatchingGraph(10));
    auto &graph = flooder.graph;
    graph.add_edge(0, 1, 10, 0);
    graph.add_edge(1, 2, 10, 0);
    graph.add_edge(2, 3, 10, 0);
    graph.add_edge(3, 4, 10, 0);
    graph.add_edge(4, 5, 10, 0);
    graph.add_edge(5, 0, 10, 0);

    auto qn = [&](int i, int t) {
        graph.nodes[i].node_event_tracker.set_desired_event(
            {&graph.nodes[i],
             cyclic_time_int{t}},
            flooder.queue);
    };

    qn(0, 10);
    qn(1, 8);
    qn(2, 5);

    GraphFillRegion gfr;
    gfr.shrink_event_tracker.set_desired_event(
        {&gfr, cyclic_time_int{70}}, flooder.queue);

    qn(4, 100);
    auto e = flooder.dequeue_valid();
    ASSERT_EQ(e.time, 5);
    ASSERT_EQ(e.tentative_event_type, LOOK_AT_NODE);
    ASSERT_EQ(e.data_look_at_node, &graph.nodes[2]);
    ASSERT_EQ(flooder.dequeue_valid().time, 8);
    ASSERT_EQ(flooder.dequeue_valid().time, 10);
    auto e6 = flooder.dequeue_valid();
    ASSERT_EQ(e6.time, 70);
    ASSERT_EQ(e6.tentative_event_type, LOOK_AT_SHRINKING_REGION);
    ASSERT_EQ(e6.data_look_at_shrinking_region, &gfr);
    ASSERT_EQ(flooder.dequeue_valid().time, 100);
}

TEST(GraphFlooder, CreateRegion) {
    GraphFlooder flooder(MatchingGraph(5));
    auto &g = flooder.graph;
    g.add_boundary_edge(0, 3, 0);
    g.add_edge(0, 1, 5, 0);
    g.add_edge(1, 2, 11, 0);
    g.add_edge(2, 3, 100, 0);
    g.add_boundary_edge(3, 1000, 0);
    flooder.create_region(&flooder.graph.nodes[0]);
    flooder.create_region(&flooder.graph.nodes[2]);
    flooder.create_region(&flooder.graph.nodes[3]);
    ASSERT_EQ(
        flooder.queue.dequeue(),
        FloodCheckEvent(&flooder.graph.nodes[0], cyclic_time_int{3}));
    ASSERT_EQ(
        flooder.queue.dequeue(),
        FloodCheckEvent(&flooder.graph.nodes[2], cyclic_time_int{11}));
    ASSERT_EQ(flooder.graph.nodes[0].region_that_arrived->shell_area[0], &flooder.graph.nodes[0]);
    ASSERT_EQ(flooder.graph.nodes[0].distance_from_source, 0);
    ASSERT_EQ(flooder.graph.nodes[0].reached_from_source, &flooder.graph.nodes[0]);
    ASSERT_EQ(
        flooder.queue.dequeue(),
        FloodCheckEvent(&flooder.graph.nodes[3], cyclic_time_int{50}));
    ASSERT_TRUE(flooder.queue.empty());
}

TEST(GraphFlooder, RegionGrowingToBoundary) {
    GraphFlooder flooder(MatchingGraph(10));
    auto &g = flooder.graph;
    g.add_boundary_edge(0, 2, 3);
    g.add_edge(0, 1, 10, 5);
    g.add_edge(1, 2, 21, 0);
    g.add_edge(2, 3, 100, 1);
    g.add_edge(3, 4, 7, 9);
    g.add_boundary_edge(4, 5, 2);
    flooder.create_region(&flooder.graph.nodes[2]);
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
    GraphFlooder flooder(MatchingGraph(10));
    auto &g = flooder.graph;
    g.add_boundary_edge(0, 2000, 3);
    g.add_edge(0, 1, 10, 5);
    g.add_edge(1, 2, 24, 0);
    g.add_edge(2, 3, 38, 1);
    g.add_edge(3, 4, 26, 9);
    g.add_boundary_edge(4, 1000, 2);
    flooder.create_region(&flooder.graph.nodes[2]);
    flooder.create_region(&flooder.graph.nodes[4]);
    auto e1 = flooder.run_until_next_mwpm_notification();
    ASSERT_EQ(flooder.queue.cur_time, 32);
    MwpmEvent e1_exp = RegionHitRegionEventData{
        flooder.graph.nodes[4].region_that_arrived,
        flooder.graph.nodes[2].region_that_arrived,
        CompressedEdge{&flooder.graph.nodes[4], &flooder.graph.nodes[2], 8}};
    ASSERT_EQ(e1, e1_exp);
}

TEST(GraphFlooder, RegionGrowingThenFrozenThenStartShrinking) {
    GraphFlooder flooder(MatchingGraph(10));
    auto &g = flooder.graph;
    g.add_boundary_edge(0, 4, 3);
    g.add_edge(0, 1, 10, 5);
    g.add_edge(1, 2, 22, 0);
    g.add_edge(2, 3, 30, 1);
    g.add_edge(3, 4, 50, 9);
    g.add_boundary_edge(4, 100, 2);
    flooder.create_region(&flooder.graph.nodes[2]);
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

    std::vector<DetectorNode *> expected_area = {
        &flooder.graph.nodes[2], &flooder.graph.nodes[1], &flooder.graph.nodes[3], &flooder.graph.nodes[0]};
    ASSERT_EQ(flooder.graph.nodes[2].region_that_arrived->shell_area, expected_area);
    flooder.set_region_shrinking(*flooder.graph.nodes[2].region_that_arrived);
    ASSERT_EQ(
        flooder.dequeue_valid(),
        FloodCheckEvent(
            flooder.graph.nodes[2].region_that_arrived,
            cyclic_time_int{80 + 4}));
}
