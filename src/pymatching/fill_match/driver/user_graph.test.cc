#include <cmath>
#include <gtest/gtest.h>

#include "pymatching/fill_match/driver/user_graph.h"

TEST(PythonAPIGraph, ConstructGraph) {
    pm::UserGraph graph;
    graph.add_or_merge_boundary_edge(0, {2}, 4.1, 0.1);
    graph.add_or_merge_boundary_edge(0, {3}, 1.0, 0.46);
    graph.add_or_merge_edge(0, 1, {0}, 2.5, 0.4);
    graph.add_or_merge_edge(0, 1, {0}, 2.1, 0.45);
    graph.add_or_merge_edge(1, 2, {1}, 3.5, 0.2);
    graph.add_or_merge_edge(2, 3, {3}, 1.8, 0.3);
    graph.add_or_merge_edge(2, 4, {4}, 2.0, 0.25);
    graph.add_or_merge_boundary_edge(2, {4}, 2.2, 0.25);
    graph.set_boundary({3, 4});
    ASSERT_EQ(graph.get_num_observables(), 5);
    ASSERT_EQ(graph.nodes.size(), 5);
    ASSERT_EQ(graph.nodes[0].neighbors[0].node, SIZE_MAX);
    ASSERT_EQ(graph.nodes[0].neighbors[0].weight, pm::merge_weights(4.1, 1.0));
    ASSERT_EQ(graph.nodes[0].neighbors[0].error_probability, 0.1 * (1 - 0.46) + 0.46 * (1 - 0.1));
    ASSERT_EQ(graph.nodes[0].neighbors[1].node, 1);
    ASSERT_EQ(graph.nodes[0].neighbors[1].weight, pm::merge_weights(2.5, 2.1));
    ASSERT_EQ(graph.nodes[0].neighbors[1].error_probability, 0.4 * (1 - 0.45) + 0.45 * (1 - 0.4));
    ASSERT_EQ(graph.nodes[1].neighbors[0].weight, pm::merge_weights(2.5, 2.1));
    ASSERT_EQ(graph.nodes[1].neighbors[0].error_probability, 0.4 * (1 - 0.45) + 0.45 * (1 - 0.4));
    ASSERT_EQ(graph.nodes[1].neighbors[0].node, 0);
    ASSERT_EQ(graph.nodes[2].neighbors[0].node, 1);
    ASSERT_EQ(graph.nodes[2].neighbors[1].node, 3);
    ASSERT_EQ(graph.nodes[0].index_of_neighbor(1), 1);
    ASSERT_EQ(graph.nodes[0].index_of_neighbor(SIZE_MAX), 0);
    ASSERT_EQ(graph.nodes[0].index_of_neighbor(3), SIZE_MAX);
    auto& mwpm = graph.get_mwpm();
    auto& g2 = mwpm.flooder.graph;
    ASSERT_EQ(g2.nodes.size(), 5);
    pm::Neighbor n = {nullptr, 4.1, {2}};
    ASSERT_EQ(g2.nodes[0].neighbors[0], nullptr);
    ASSERT_EQ(g2.nodes[0].neighbor_weights[0], (pm::weight_int)(pm::merge_weights(4.1, 1.0) * mwpm.flooder.graph.normalising_constant));
    ASSERT_EQ(g2.nodes[0].neighbor_observables[0], 1 << 2);
    ASSERT_EQ(g2.nodes[0].neighbors[1], &g2.nodes[1]);
    ASSERT_EQ(g2.nodes[0].neighbor_weights[1], 2 * (pm::weight_int)(pm::merge_weights(2.5, 2.1) * mwpm.flooder.graph.normalising_constant / 2));
    ASSERT_EQ(g2.nodes[0].neighbor_observables[1], 1 << 0);
    ASSERT_EQ(g2.nodes[2].neighbors[0], nullptr);
    ASSERT_EQ(g2.nodes[2].neighbor_weights[0], 2 * (pm::weight_int)(1.8 * mwpm.flooder.graph.normalising_constant / 2));
    ASSERT_EQ(g2.nodes[2].neighbor_observables[0], 1 << 3);
    ASSERT_EQ(g2.nodes[2].neighbors[1], &g2.nodes[1]);
    ASSERT_EQ(g2.nodes[2].neighbor_weights[1], 2 * (pm::weight_int)(3.5 * mwpm.flooder.graph.normalising_constant / 2));
    ASSERT_EQ(g2.nodes[2].neighbor_observables[1], 1 << 1);
    ASSERT_EQ(g2.nodes[2].neighbors.size(), 2);
    ASSERT_EQ(g2.nodes[3].neighbors.size(), 0);
}

TEST(PythonApiGraph, AddNoise) {
    pm::UserGraph graph;
    graph.add_or_merge_boundary_edge(0, {0}, 1, 1);
    graph.add_or_merge_edge(0, 1, {1}, 1, 0);
    graph.add_or_merge_edge(1, 2, {2}, 1, 0);
    graph.add_or_merge_edge(2, 3, {3}, 1, 1);
    graph.add_or_merge_edge(3, 4, {4}, 1, 0);
    graph.add_or_merge_edge(4, 5, {5}, 1, 0);
    graph.add_or_merge_edge(5, 6, {6}, 1, 1);
    graph.add_or_merge_edge(6, 7, {7}, 1, 1);
    graph.set_boundary({7});
    std::vector<uint8_t> observables(graph.get_num_observables());
    std::vector<uint8_t> syndrome(graph.get_num_nodes());
    graph.add_noise(observables.data(), syndrome.data());
    std::vector<uint8_t> expected_observables = {1, 0, 0, 1, 0, 0, 1, 1};
    ASSERT_EQ(observables, expected_observables);
    std::vector<uint8_t> expected_syndrome = {1, 0, 1, 1, 0, 1, 0, 0};
    ASSERT_EQ(syndrome, expected_syndrome);
}
