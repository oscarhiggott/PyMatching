#include "pymatching/fill_match/driver/python_api_graph.h"

#include <cmath>
#include <gtest/gtest.h>

TEST(PythonAPIGraph, ConstructGraph) {
    pm::UserGraph graph;
    graph.add_boundary_edge(0, {2}, 4.1, 0.1);
    graph.add_edge(0, 1, {0}, 2.5, 0.4);
    graph.add_edge(1, 2, {1}, 3.5, 0.2);
    graph.add_edge(2, 3, {3}, 1.8, 0.3);
    graph.set_boundary({3});
    ASSERT_EQ(graph.get_num_observables(), 4);
    ASSERT_EQ(graph.nodes.size(), 4);
    ASSERT_EQ(graph.nodes[0].neighbors[0].node, SIZE_MAX);
    ASSERT_EQ(graph.nodes[0].neighbors[1].node, 1);
    ASSERT_EQ(graph.nodes[1].neighbors[0].node, 0);
    ASSERT_EQ(graph.nodes[2].neighbors[0].node, 1);
    ASSERT_EQ(graph.nodes[2].neighbors[1].node, 3);
    auto g2 = graph.to_intermediate_weighted_graph();
    ASSERT_EQ(g2.nodes.size(), 4);
    pm::Neighbor n = {nullptr, 4.1, {2}};
    ASSERT_EQ(g2.nodes[0][0], n);
    pm::Neighbor n2 = {&g2.nodes[1], 2.5, {0}};
    ASSERT_EQ(g2.nodes[0][1], n2);
    pm::Neighbor n3 = {nullptr, 1.8, {3}};
    ASSERT_EQ(g2.nodes[2][1], n3);
    pm::Neighbor n4 = {&g2.nodes[1], 3.5, {1}};
    ASSERT_EQ(g2.nodes[2][0].node, &g2.nodes[1]);
    ASSERT_EQ(g2.nodes[3].size(), 0);
}
