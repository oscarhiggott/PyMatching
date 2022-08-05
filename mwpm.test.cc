#include <gtest/gtest.h>
#include "mwpm.h"

TEST(Mwpm, BlossomCreatedThenShattered) {
    auto g = pm::Graph(10);
    g.add_edge(0, 1, 10, 1);
    g.add_edge(1, 4, 20, 2);
    g.add_edge(4, 3, 20, 3);
    g.add_edge(3, 2, 12, 4);
    g.add_edge(0, 2, 16, 5);
    g.add_edge(4, 5, 50, 6);
    g.add_edge(2, 6, 100, 7);
    g.add_boundary_edge(5, 36, 8);
    auto flooder = pm::GraphFlooder(g);
    auto mwpm = pm::Mwpm(flooder);
    for (size_t i = 0; i < 7; i++)
        mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[i]);
    auto& ns = mwpm.flooder.graph.nodes;
    auto e1 = mwpm.flooder.next_event();
    mwpm.process_event(e1);
    auto e2 = mwpm.flooder.next_event();
    mwpm.process_event(e2);
    auto e3 = mwpm.flooder.next_event();
    mwpm.process_event(e3);
    auto e4 = mwpm.flooder.next_event();
    mwpm.process_event(e4);
    auto e5 = mwpm.flooder.next_event();
    ASSERT_EQ(e5, pm::MwpmEvent(
            mwpm.flooder.graph.nodes[0].region_that_arrived,
            mwpm.flooder.graph.nodes[2].region_that_arrived,
            {&mwpm.flooder.graph.nodes[0], &mwpm.flooder.graph.nodes[2], 5}
            ));
    mwpm.process_event(e5);
    // Region 5 matches with blossom containing {0, 1, 4, 3, 2}
    auto e6 = mwpm.flooder.next_event();
    mwpm.process_event(e6);
    auto blossom = mwpm.flooder.graph.nodes[5].region_that_arrived->match.region;
    ASSERT_EQ(blossom->blossom_children.size(), 5);
    std::vector<pm::RegionEdge> expected_blossom_children = {
            {ns[2].region_that_arrived, {&ns[2], &ns[3], 4}},
            {ns[3].region_that_arrived, {&ns[3], &ns[4], 3}},
            {ns[4].region_that_arrived, {&ns[4], &ns[1], 2}},
            {ns[1].region_that_arrived, {&ns[1], &ns[0], 1}},
            {ns[0].region_that_arrived, {&ns[0], &ns[2], 5}},
    };
    ASSERT_EQ(blossom->blossom_children, expected_blossom_children);
    ASSERT_EQ(blossom->radius, pm::Varying(8 << 2));
    ASSERT_EQ(mwpm.flooder.time, 25);
    // Region 6 collides with matched blossom
    auto e7 = mwpm.flooder.next_event();
    mwpm.process_event(e7);
    ASSERT_EQ(mwpm.flooder.time, 83);
    ASSERT_EQ(
            mwpm.flooder.graph.nodes[6].region_that_arrived->alt_tree_node->children[0].alt_tree_node->inner_region,
            blossom
              );
    // Blossom shatters
    auto e8 = mwpm.flooder.next_event();
    ASSERT_EQ(mwpm.flooder.time, 91);
    ASSERT_EQ(e8.event_type, pm::BLOSSOM_SHATTER);
    mwpm.process_event(e8);
    ASSERT_EQ(ns[0].region_that_arrived->match.region,
              ns[1].region_that_arrived);
    ASSERT_EQ(ns[1].region_that_arrived->match.region,
              ns[0].region_that_arrived);
    auto e9 = mwpm.flooder.next_event();
    ASSERT_EQ(mwpm.flooder.time, 94);
    ASSERT_EQ(e9.event_type, pm::REGION_HIT_BOUNDARY);
    mwpm.process_event(e9);
    auto e10 = mwpm.flooder.next_event();
    ASSERT_EQ(e10.event_type, pm::NO_EVENT);
    ASSERT_EQ(ns[3].region_that_arrived->blossom_parent, nullptr);
    auto regions_matched = [&ns](size_t i, size_t j) {
        bool regions_correct = ns[i].region_that_arrived->match.region == ns[j].region_that_arrived &&
                                ns[j].region_that_arrived->match.region == ns[i].region_that_arrived;
        bool edge_locs_correct = ns[i].region_that_arrived->match.edge.loc_from == &ns[i] &&
                                  ns[i].region_that_arrived->match.edge.loc_to == &ns[j] &&
                                  ns[j].region_that_arrived->match.edge.loc_from == &ns[j] &&
                                  ns[j].region_that_arrived->match.edge.loc_to == &ns[i];
        return regions_correct && edge_locs_correct;
    };
    ASSERT_TRUE(regions_matched(4, 3));
    ASSERT_TRUE(regions_matched(2, 6));
    ASSERT_TRUE(regions_matched(1, 0));
    ASSERT_EQ(ns[5].region_that_arrived->match, pm::Match(
            nullptr, {&ns[5], nullptr, 8}
            ));
}
