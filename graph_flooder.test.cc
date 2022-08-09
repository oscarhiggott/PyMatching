#include "gtest/gtest.h"
#include "graph_flooder.h"
#include "events.h"
#include "mwpm.h"

TEST(GraphFlooder, PriorityQueue){
    pm::MatchingGraph graph(10);
    pm::GraphFillRegion gfr;
    pm::GraphFlooder flooder(graph);
    auto e1 = new pm::TentativeEvent(&graph.nodes[0], 0, &graph.nodes[1], 1, 10);
    flooder.queue.push(e1);
    auto e2 = new pm::TentativeEvent(&graph.nodes[1], 0, &graph.nodes[2], 1, 8);
    flooder.queue.push(e2);
    auto e3 = new pm::TentativeEvent(&graph.nodes[2], 0, &graph.nodes[3], 1, 5);
    flooder.queue.push(e3);
    auto e4 = new pm::TentativeEvent(&gfr, 70);
    flooder.queue.push(e4);
    auto e5 = new pm::TentativeEvent(&graph.nodes[4], 0, &graph.nodes[5], 1, 100);
    flooder.queue.push(e5);
    auto e = flooder.queue.top();
    ASSERT_EQ(e->time, 5);
    ASSERT_EQ(e->tentative_event_type, pm::INTERACTION);
    ASSERT_EQ(e->neighbor_interaction_event_data.detector_node_1, &graph.nodes[2]);
    flooder.queue.pop();
    ASSERT_EQ(flooder.queue.top()->time, 8);
    flooder.queue.pop();
    ASSERT_EQ(flooder.queue.top()->time, 10);
    flooder.queue.pop();
    auto e6 = flooder.queue.top();
    ASSERT_EQ(e6->time, 70);
    ASSERT_EQ(e6->tentative_event_type, pm::SHRINKING);
    ASSERT_EQ(e6->region_shrink_event_data.region, &gfr);
    ASSERT_EQ(e6->is_invalidated, false);
    flooder.queue.pop();
    ASSERT_EQ(flooder.queue.top()->time, 100);
}


TEST(GraphFlooder, CreateRegion) {
    auto g = pm::MatchingGraph(5);
    g.add_boundary_edge(0, 3, 0);
    g.add_edge(0, 1, 5, 0);
    g.add_edge(1, 2, 11, 0);
    g.add_edge(2, 3, 100, 0);
    g.add_boundary_edge(3, 1000, 0);
    pm::GraphFlooder flooder(g);
    flooder.create_region(&flooder.graph.nodes[0]);
    flooder.create_region(&flooder.graph.nodes[2]);
    flooder.create_region(&flooder.graph.nodes[3]);
    auto e1 = flooder.queue.top();
    ASSERT_EQ(*e1, pm::TentativeEvent(
            &flooder.graph.nodes[0], 0,
            nullptr, -1, 3
            ));
    flooder.queue.pop();
    auto e2 = flooder.queue.top();
    ASSERT_EQ(*e2,
            pm::TentativeEvent(
                    &flooder.graph.nodes[0], 1,
                    &flooder.graph.nodes[1], 0, 5
                    )
            );
    flooder.queue.pop();
    ASSERT_EQ(flooder.graph.nodes[0].region_that_arrived->shell_area[0], &flooder.graph.nodes[0]);
    ASSERT_EQ(flooder.graph.nodes[0].distance_from_source, 0);
    ASSERT_EQ(flooder.graph.nodes[0].reached_from_source, &flooder.graph.nodes[0]);
    auto e3 = flooder.queue.top();
    ASSERT_EQ(*e3,
              pm::TentativeEvent(
                      &flooder.graph.nodes[2], 0,
                      &flooder.graph.nodes[1], 1,
                      11
                      )
              );
    flooder.queue.pop();
    auto e4 = flooder.queue.top();
    ASSERT_EQ(*e4,
              pm::TentativeEvent(
                      &flooder.graph.nodes[3], 0,
                      &flooder.graph.nodes[2], 1, 50
                      )
              );
    flooder.queue.pop();
    auto e5 = flooder.queue.top();
    auto e5_exp = pm::TentativeEvent(
            &flooder.graph.nodes[2], 1,
            &flooder.graph.nodes[3], 0,
            100
    );
    e5_exp.is_invalidated = true;
    ASSERT_EQ(*e5, e5_exp);
    flooder.queue.pop();
    auto e6 = flooder.queue.top();
    ASSERT_EQ(*e6,
              pm::TentativeEvent(
                      &flooder.graph.nodes[3], 1,
                      nullptr, -1, 1000
                      )
              );
    ASSERT_FALSE(flooder.queue.empty());
    flooder.queue.pop();
    ASSERT_TRUE(flooder.queue.empty());
}


TEST(GraphFlooder, RegionGrowingToBoundary){
    auto g = pm::MatchingGraph(10);
    g.add_boundary_edge(0, 2, 3);
    g.add_edge(0, 1, 10, 5);
    g.add_edge(1, 2, 21, 0);
    g.add_edge(2, 3, 100, 1);
    g.add_edge(3, 4, 7, 9);
    g.add_boundary_edge(4, 5, 2);
    pm::GraphFlooder flooder(g);
    flooder.create_region(&flooder.graph.nodes[2]);
    auto e1 = flooder.next_event();
    auto e1_exp = pm::MwpmEvent(
            flooder.graph.nodes[2].region_that_arrived,
            pm::CompressedEdge(&flooder.graph.nodes[2], nullptr, 6)
    );
    ASSERT_EQ(e1, e1_exp);
    ASSERT_EQ(*flooder.queue.top(),
              pm::TentativeEvent(
                      &flooder.graph.nodes[2],
                      1,
                      &flooder.graph.nodes[3],
                      0,
                      100
                      )
              );
    auto e2 = flooder.next_event();
    auto e2_exp = pm::MwpmEvent(
            flooder.graph. nodes[2].region_that_arrived,
            pm::CompressedEdge(&flooder.graph.nodes[2], nullptr, 10)
    );
    ASSERT_EQ(
            e2,
            e2_exp
            );
    auto e3 = flooder.next_event();
    ASSERT_EQ(e3.event_type, pm::NO_EVENT);
}


TEST(GraphFlooder, RegionHitRegion) {
    auto g = pm::MatchingGraph(10);
    g.add_boundary_edge(0, 2000, 3);
    g.add_edge(0, 1, 10, 5);
    g.add_edge(1, 2, 24, 0);
    g.add_edge(2, 3, 38, 1);
    g.add_edge(3, 4, 26, 9);
    g.add_boundary_edge(4, 1000, 2);
    pm::GraphFlooder flooder(g);
    flooder.create_region(&flooder.graph.nodes[2]);
    flooder.create_region(&flooder.graph.nodes[4]);
    auto e1 = flooder.next_event();
    ASSERT_EQ(flooder.time, 32);
    auto e1_exp = pm::MwpmEvent(
            flooder.graph.nodes[4].region_that_arrived,
            flooder.graph.nodes[2].region_that_arrived,
            pm::CompressedEdge(
                    &flooder.graph.nodes[4],
                    &flooder.graph.nodes[2],
                    8
                    )
            );
    ASSERT_EQ(e1, e1_exp);
}


TEST(GraphFlooder, RegionGrowingThenFrozenThenStartShrinking) {
    auto g = pm::MatchingGraph(10);
    g.add_boundary_edge(0, 4, 3);
    g.add_edge(0, 1, 10, 5);
    g.add_edge(1, 2, 22, 0);
    g.add_edge(2, 3, 30, 1);
    g.add_edge(3, 4, 50, 9);
    g.add_boundary_edge(4, 100, 2);
    pm::GraphFlooder flooder(g);
    flooder.create_region(&flooder.graph.nodes[2]);
    auto e1 = flooder.next_event();
    ASSERT_EQ(
            e1,
            pm::MwpmEvent(
                    flooder.graph.nodes[2].region_that_arrived,
                    pm::CompressedEdge(
                            &flooder.graph.nodes[2],
                            nullptr,
                            6
                            )
                    )
              );
    ASSERT_EQ(flooder.time, 36);
    flooder.set_region_frozen(*flooder.graph.nodes[2].region_that_arrived);
    auto e2 = flooder.next_event();
    ASSERT_EQ(flooder.time, 36);
    ASSERT_EQ(e2.event_type, pm::NO_EVENT);
    ASSERT_TRUE(flooder.queue.empty());
    std::vector<pm::DetectorNode*> expected_area = {&flooder.graph.nodes[2], &flooder.graph.nodes[1],
                                                    &flooder.graph.nodes[3], &flooder.graph.nodes[0]};
    ASSERT_EQ(flooder.graph.nodes[2].region_that_arrived->shell_area, expected_area);
    flooder.set_region_shrinking(*flooder.graph.nodes[2].region_that_arrived);
    ASSERT_EQ(
            *flooder.queue.top(),
            pm::TentativeEvent(
                    flooder.graph.nodes[2].region_that_arrived,
                    40
                    )
              );
}


TEST(GraphFlooder, TwoRegionsGrowingThenMatching) {
    auto g = pm::MatchingGraph(10);
    g.add_boundary_edge(0, 4, 3);
    g.add_edge(0, 1, 100, 5);
    g.add_edge(1, 2, 22, 5);
    g.add_edge(2, 3, 30, 1);
    g.add_edge(3, 4, 10, 9);
    g.add_boundary_edge(4, 1000, 2);
    pm::GraphFlooder flooder(g);
    pm::Mwpm mwpm(flooder);
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[1]);
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[3]);
    auto e1 = mwpm.flooder.next_event();
    auto e1_expected = pm::MwpmEvent(
            mwpm.flooder.graph.nodes[1].region_that_arrived,
            mwpm.flooder.graph.nodes[3].region_that_arrived,
            pm::CompressedEdge(
                    &mwpm.flooder.graph.nodes[1],
                    &mwpm.flooder.graph.nodes[3],
                    4
                    )
            );
    ASSERT_EQ(e1, e1_expected);
    mwpm.process_event(e1);
    auto e2 = mwpm.flooder.next_event();
    ASSERT_EQ(e2.event_type, pm::NO_EVENT);
    ASSERT_TRUE(mwpm.flooder.queue.empty());
    std::vector<pm::DetectorNode*> expected_area_1 = {
            &mwpm.flooder.graph.nodes[1], &mwpm.flooder.graph.nodes[2]
    };
    ASSERT_EQ(
            mwpm.flooder.graph.nodes[1].region_that_arrived->shell_area,
            expected_area_1
            );
    std::vector<pm::DetectorNode*> expected_area_3 = {
            &mwpm.flooder.graph.nodes[3], &mwpm.flooder.graph.nodes[4]
    };
    ASSERT_EQ(
            mwpm.flooder.graph.nodes[3].region_that_arrived->shell_area,
            expected_area_3
    );
    ASSERT_EQ(mwpm.flooder.graph.nodes[1].region_that_arrived->radius, pm::Varying(26 << 2));
    ASSERT_EQ(mwpm.flooder.graph.nodes[3].region_that_arrived->radius, pm::Varying(26 << 2));
}

TEST(GraphFlooder, RegionHittingMatchThenMatchedToOtherRegion) {
    auto g = pm::MatchingGraph(10);
    g.add_boundary_edge(0, 1000, 3);
    g.add_edge(0, 1, 8, 5);
    g.add_edge(1, 2, 10, 5);
    g.add_edge(2, 3, 2, 1);
    g.add_edge(3, 4, 4, 9);
    g.add_edge(4, 5, 20, 2);
    g.add_edge(5, 6, 36, 3);
    g.add_boundary_edge(6, 1000, 2);
    pm::GraphFlooder flooder(g);
    pm::Mwpm mwpm(flooder);
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[1]);
    auto r1 = mwpm.flooder.graph.nodes[1].region_that_arrived;
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[4]);
    auto r4 = mwpm.flooder.graph.nodes[4].region_that_arrived;
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[5]);
    auto r5 = mwpm.flooder.graph.nodes[5].region_that_arrived;
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[6]);
    auto r6 = mwpm.flooder.graph.nodes[6].region_that_arrived;
    auto e1 = mwpm.flooder.next_event();
    auto e1_expected = pm::MwpmEvent(
            r4,
            r1,
            pm::CompressedEdge(
                    &mwpm.flooder.graph.nodes[4],
                    &mwpm.flooder.graph.nodes[1],
                    13
                    )
            );
    ASSERT_EQ(e1, e1_expected);
    mwpm.process_event(e1);
    auto e2 = mwpm.flooder.next_event();
    auto e2_expected = pm::MwpmEvent(
            r4,
            r5,
            pm::CompressedEdge(
                    &mwpm.flooder.graph.nodes[4],
                    &mwpm.flooder.graph.nodes[5],
                    2
                    )
            );
    ASSERT_EQ(e2, e2_expected);
    ASSERT_EQ(mwpm.flooder.time, 12);
    mwpm.process_event(e2);
    auto e3 = mwpm.flooder.next_event();
    auto e3_expected = pm::MwpmEvent(
            r6,
            r5,
            pm::CompressedEdge(
                    &mwpm.flooder.graph.nodes[6],
                    &mwpm.flooder.graph.nodes[5],
                    3
                    )
            );
    ASSERT_EQ(
        e3,
        e3_expected
            );
    ASSERT_EQ(mwpm.flooder.time, 18);
    mwpm.process_event(e3);
    auto e4 = mwpm.flooder.next_event();
    ASSERT_EQ(e4.event_type, pm::NO_EVENT);
    ASSERT_EQ(r1->radius, pm::Varying(14 << 2));
    std::vector<pm::DetectorNode*> area_1 = {
            &mwpm.flooder.graph.nodes[1],
            &mwpm.flooder.graph.nodes[0],
            &mwpm.flooder.graph.nodes[2],
            &mwpm.flooder.graph.nodes[3]
    };
    ASSERT_EQ(r1->shell_area, area_1);
    ASSERT_EQ(r4->radius, pm::Varying(2 << 2));
    std::vector<pm::DetectorNode*> area_4 = {&mwpm.flooder.graph.nodes[4]};
    ASSERT_EQ(r4->shell_area, area_4);
    for (auto ri : {r1, r4, r5, r6})
        ASSERT_EQ(ri->alt_tree_node, nullptr);
    ASSERT_EQ(
            r1->match,
            pm::Match(
                    r4,
                    pm::CompressedEdge(
                            &mwpm.flooder.graph.nodes[1],
                            &mwpm.flooder.graph.nodes[4],
                            13
                            )
                    )
              );
    ASSERT_EQ(
            r4->match,
            pm::Match(
                    r1,
                    pm::CompressedEdge(
                            &mwpm.flooder.graph.nodes[4],
                            &mwpm.flooder.graph.nodes[1],
                            13
                    )
            )
    );
    ASSERT_EQ(
            r5->match,
            pm::Match(
                    r6,
                    pm::CompressedEdge(
                            &mwpm.flooder.graph.nodes[5],
                            &mwpm.flooder.graph.nodes[6],
                            3
                    )
            )
    );
    ASSERT_EQ(
            r6->match,
            pm::Match(
                    r5,
                    pm::CompressedEdge(
                            &mwpm.flooder.graph.nodes[6],
                            &mwpm.flooder.graph.nodes[5],
                            3
                    )
            )
    );
}

TEST(GraphFlooder, RegionHittingMatchFormingBlossomThenMatchingToBoundary) {
    size_t num_nodes = 100;
    auto g = pm::MatchingGraph(num_nodes);
    g.add_boundary_edge(0, 2, 1);
    for (size_t i = 0; i < num_nodes - 1; i++)
        g.add_edge(i, i + 1, 2, i);
    g.add_boundary_edge(num_nodes - 1, 2, 2);
    pm::GraphFlooder flooder(g);
    pm::Mwpm mwpm(flooder);
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[40]);
    auto r40 = mwpm.flooder.graph.nodes[40].region_that_arrived;
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[42]);
    auto r42 = mwpm.flooder.graph.nodes[42].region_that_arrived;
    mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[50]);
    auto r50 = mwpm.flooder.graph.nodes[50].region_that_arrived;
    auto e1 = mwpm.flooder.next_event();
    mwpm.process_event(e1);
    auto e2 = mwpm.flooder.next_event();
    mwpm.process_event(e2);
    auto e3 = mwpm.flooder.next_event();
    mwpm.process_event(e3);
    auto e4 = mwpm.flooder.next_event();
    mwpm.process_event(e4);
    ASSERT_EQ(mwpm.flooder.time, 94);
    std::vector<pm::RegionEdge> expected_blossom_regions = {
            {r40,
                {&mwpm.flooder.graph.nodes[40], &mwpm.flooder.graph.nodes[42], 40 ^ 41 }
                },
            {
                    r42,
                    {
                            &mwpm.flooder.graph.nodes[42], &mwpm.flooder.graph.nodes[50],
                            42 ^ 43 ^ 44 ^ 45 ^ 46 ^ 47 ^ 48 ^ 49
                    }
            },{
                r50,
                {&mwpm.flooder.graph.nodes[50], &mwpm.flooder.graph.nodes[40],
                    40 ^ 41 ^ 42 ^ 43 ^ 44 ^ 45 ^ 46 ^ 47 ^ 48 ^ 49
                    }},

    };
    auto blossom = mwpm.flooder.graph.nodes[40].region_that_arrived->blossom_parent;
    auto actual_blossom_regions =blossom->blossom_children;
    ASSERT_EQ(actual_blossom_regions, expected_blossom_regions);
    pm::obs_int obs = 0;
    for (int i = 0; i < 40; i++)
        obs ^= i;
    ASSERT_EQ(
            blossom->match,
            pm::Match(
                    nullptr, {&mwpm.flooder.graph.nodes[40], nullptr, obs ^ 1}
                    )
              );
}
