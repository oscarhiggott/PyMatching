#include "pymatching/fill_match/flooder/graph.h"

#include <gtest/gtest.h>

#include "pymatching/fill_match/flooder/graph_fill_region.test.h"

TEST(Graph, AddEdge) {
    pm::MatchingGraph g(4, 64);
    g.add_edge(0, 1, 2, {0});
    g.add_edge(1, 2, 3, {0, 2});
    g.add_edge(0, 3, 10, {1, 3});
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
    pm::MatchingGraph g(6, 64);
    g.add_edge(0, 1, 2, {0, 1});
    g.add_boundary_edge(0, 7, {2});
    g.add_boundary_edge(5, 10, {0, 1, 3});
    ASSERT_EQ(g.nodes[0].neighbors[0], nullptr);
    ASSERT_EQ(g.nodes[0].neighbors[1], &g.nodes[1]);
    ASSERT_EQ(g.nodes[5].neighbors[0], nullptr);
}
