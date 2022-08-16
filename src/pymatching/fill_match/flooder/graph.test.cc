#include "pymatching/fill_match/flooder/graph.h"

#include <gtest/gtest.h>

#include "pymatching/fill_match/flooder/graph_fill_region.test.h"

TEST(Graph, AddEdge) {
    pm::MatchingGraph g(4);
    g.add_edge(0, 1, 2, 1);
    g.add_edge(1, 2, 3, 5);
    g.add_edge(0, 3, 10, 10);
    ASSERT_EQ(g.nodes[0].neighbors[0], &g.nodes[1]);
    ASSERT_EQ(g.nodes[1].neighbors[0], &g.nodes[0]);
    ASSERT_EQ(g.nodes[3].neighbors[0], &g.nodes[0]);
    ASSERT_EQ(g.nodes[0].neighbors[1], &g.nodes[3]);
    ASSERT_EQ(g.nodes[0].neighbor_weights[0], 2);
    ASSERT_EQ(g.nodes[0].neighbor_weights[1], 10);
    ASSERT_EQ(g.nodes[1].neighbor_observables[0], 1);
    ASSERT_EQ(g.nodes[3].neighbor_observables[0], 10);
}

TEST(Graph, AddBoundaryEdge) {
    pm::MatchingGraph g(6);
    g.add_edge(0, 1, 2, 3);
    g.add_boundary_edge(0, 7, 4);
    g.add_boundary_edge(5, 10, 11);
    ASSERT_EQ(g.nodes[0].neighbors[0], nullptr);
    ASSERT_EQ(g.nodes[0].neighbors[1], &g.nodes[1]);
    ASSERT_EQ(g.nodes[5].neighbors[0], nullptr);
}

TEST(Graph, TotalRadius) {
    pm::GraphFillTestData d(10);
    auto x = d.b(-1,
                 -1,
                 {d.b(0, 1, {d.b(3, 4, {}, false), d.b(4, 5, {}, false), d.b(5, 3, {}, false)}, false),
                  d.b(1, 2, {}, false),
                  d.b(2, 0, {}, false)},
                 true)
                 .region;
    d.detectors[3].reached_from_source = &d.detectors[3];
    d.detectors[6].reached_from_source = &d.detectors[3];
    d.detectors[3].region_that_arrived = x->blossom_children[0].region->blossom_children[0].region;
    d.detectors[3].region_that_arrived_top = d.detectors[3].region_that_arrived->blossom_parent_top;
    d.detectors[6].region_that_arrived = x;
    d.detectors[6].region_that_arrived_top = d.detectors[6].region_that_arrived->blossom_parent_top;
    x->blossom_children[0].region->blossom_children[0].region->radius = pm::Varying32(5 << 2);
    x->blossom_children[0].region->radius = pm::Varying32(6 << 2);
    x->radius = pm::Varying32((1 << 2) + 1);
    d.detectors[6].distance_from_source = 20;
    ASSERT_EQ(d.detectors[6].total_radius(), pm::Varying32((12 << 2) + 1));
    ASSERT_EQ(d.detectors[6].local_radius(), pm::Varying32((-8 << 2) + 1));
    ASSERT_EQ(d.detectors[6].local_radius().get_distance_at_time(1000), 992);
}
