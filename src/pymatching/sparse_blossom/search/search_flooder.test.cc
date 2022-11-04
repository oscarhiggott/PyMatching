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

#include "pymatching/sparse_blossom/search/search_flooder.h"

#include "gtest/gtest.h"

TEST(SearchFlooder, RepCodeDetectorSearch) {
    size_t num_nodes = 40;
    auto flooder = pm::SearchFlooder(pm::SearchGraph(num_nodes));
    auto& g = flooder.graph;
    g.add_boundary_edge(0, 2, {0});
    for (size_t i = 0; i < num_nodes - 1; i++)
        g.add_edge(i, i + 1, 2, {i + 1});

    auto collision_edge = flooder.run_until_collision(&g.nodes[1], &g.nodes[20]);
    ASSERT_EQ(collision_edge.detector_node, &g.nodes[10]);
    ASSERT_EQ(collision_edge.neighbor_index, 1);
    ASSERT_EQ(flooder.target_type, pm::DETECTOR_NODE);
    std::vector<uint8_t> observables(num_nodes, 0);
    pm::total_weight_int weight = 0;
    flooder.iter_edges_tracing_back_from_collision_edge(collision_edge, [&](const pm::SearchGraphEdge& e) {
        auto& obs = e.detector_node->neighbor_observable_indices[e.neighbor_index];
        for (auto i : obs)
            *(observables.data() + i) ^= 1;
        weight += e.detector_node->neighbor_weights[e.neighbor_index];
    });
    std::vector<uint8_t> expected_obs(num_nodes, 0);
    for (size_t i = 1; i < 20; i++)
        expected_obs[i + 1] ^= 1;
    ASSERT_EQ(observables, expected_obs);
    ASSERT_EQ(weight, 38);
    flooder.reset_graph();
    for (auto& n : g.nodes) {
        ASSERT_EQ(n.index_of_predecessor, SIZE_MAX);
        ASSERT_EQ(n.reached_from_source, nullptr);
        ASSERT_EQ(n.node_event_tracker.desired_time, pm::cyclic_time_int{0});
        ASSERT_EQ(n.node_event_tracker.queued_time, pm::cyclic_time_int{0});
        ASSERT_EQ(n.node_event_tracker.has_desired_time, false);
        ASSERT_EQ(n.node_event_tracker.has_queued_time, false);
    }
    ASSERT_EQ(flooder.reached_nodes.size(), 0);
    flooder.reset();
}

TEST(SearchFlooder, RepCodeBoundarySearch) {
    size_t num_nodes = 20;
    auto flooder = pm::SearchFlooder(pm::SearchGraph(num_nodes));
    auto& g = flooder.graph;
    g.add_boundary_edge(0, 2, {0});
    for (size_t i = 0; i < num_nodes - 1; i++)
        g.add_edge(i, i + 1, 2, {i + 1});

    auto collision_edge = flooder.run_until_collision(&g.nodes[6], nullptr);
    ASSERT_EQ(collision_edge.detector_node, &g.nodes[0]);
    ASSERT_EQ(collision_edge.neighbor_index, 0);
    ASSERT_EQ(flooder.target_type, pm::BOUNDARY);
    std::vector<uint8_t> observables(num_nodes, 0);
    pm::total_weight_int weight = 0;
    flooder.iter_edges_tracing_back_from_collision_edge(collision_edge, [&](const pm::SearchGraphEdge& e) {
        auto& obs = e.detector_node->neighbor_observable_indices[e.neighbor_index];
        for (auto i : obs)
            *(observables.data() + i) ^= 1;
        weight += e.detector_node->neighbor_weights[e.neighbor_index];
    });
    std::vector<uint8_t> expected_obs(num_nodes, 0);
    for (size_t i = 0; i < 7; i++)
        expected_obs[i] ^= 1;
    ASSERT_EQ(observables, expected_obs);
    ASSERT_EQ(weight, 14);
    flooder.reset();
}

TEST(SearchFlooder, RepCodeSourceToDestPath) {
    size_t num_nodes = 20;
    auto flooder = pm::SearchFlooder(pm::SearchGraph(num_nodes));
    auto& g = flooder.graph;
    g.add_boundary_edge(0, 2, {0});
    for (size_t i = 0; i < num_nodes - 1; i++)
        g.add_edge(i, i + 1, 2, {i + 1});
    size_t i_start = 10;
    for (size_t i_end = 14; i_end < 18; i_end++) {
        std::vector<pm::SearchGraphEdge> edges;
        flooder.iter_edges_on_shortest_path_from_source(i_start, i_end, [&](const pm::SearchGraphEdge& edge) {
            edges.push_back(edge);
        });

        std::vector<size_t> node_indices;
        for (size_t i = 0; i < edges.size(); i++) {
            node_indices.push_back(edges[i].detector_node - &flooder.graph.nodes[0]);
            if (i + 1 < edges.size())
                ASSERT_EQ(edges[i].detector_node->neighbors[edges[i].neighbor_index], edges[i + 1].detector_node);
        }
        std::vector<size_t> expected_node_indices;
        for (size_t i = i_start; i < i_end; i++)
            expected_node_indices.push_back(i);
        ASSERT_EQ(node_indices, expected_node_indices);
    }

    std::vector<pm::SearchGraphEdge> edges;
    flooder.iter_edges_on_shortest_path_from_source(4, SIZE_MAX, [&](const pm::SearchGraphEdge& edge) {
        edges.push_back(edge);
    });
    std::vector<size_t> node_indices;
    for (size_t i = 0; i < edges.size(); i++)
        node_indices.push_back(edges[i].detector_node - &flooder.graph.nodes[0]);

    std::vector<size_t> expected_node_indices = {4, 3, 2, 1, 0};
    ASSERT_EQ(node_indices, expected_node_indices);
}
