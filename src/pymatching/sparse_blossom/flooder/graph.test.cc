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

#include "pymatching/sparse_blossom/flooder/graph.h"

#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/flooder/graph_fill_region.test.h"

TEST(Graph, AddEdge) {
    pm::MatchingGraph g(4, 64);
    g.add_edge(0, 1, -2, {0}, {});
    g.add_edge(1, 2, -3, {0, 2}, {});
    g.add_edge(0, 3, 10, {1, 3}, {});
    ASSERT_EQ(g.nodes[0].neighbors[0], &g.nodes[1]);
    ASSERT_EQ(g.nodes[1].neighbors[0], &g.nodes[0]);
    ASSERT_EQ(g.nodes[3].neighbors[0], &g.nodes[0]);
    ASSERT_EQ(g.nodes[0].neighbors[1], &g.nodes[3]);
    ASSERT_EQ(g.nodes[0].neighbor_weights[0], 2);
    ASSERT_EQ(g.nodes[0].neighbor_weights[1], 10);
    ASSERT_EQ(g.nodes[1].neighbor_observables[0], 1);
    ASSERT_EQ(g.nodes[3].neighbor_observables[0], 10);
    std::set<size_t> expected_obs = {2};
    ASSERT_EQ(g.negative_weight_observables_set, expected_obs);
    std::set<size_t> expected_detection_events = {0, 2};
    ASSERT_EQ(g.negative_weight_detection_events_set, expected_detection_events);
}

TEST(Graph, AddBoundaryEdge) {
    pm::MatchingGraph g(6, 64);
    g.add_edge(0, 1, 2, {0, 1}, {});
    g.add_boundary_edge(0, -7, {2}, {});
    g.add_boundary_edge(5, 10, {0, 1, 3}, {});
    ASSERT_EQ(g.nodes[0].neighbors[0], nullptr);
    ASSERT_EQ(g.nodes[0].neighbor_weights[0], 7);
    ASSERT_EQ(g.nodes[0].neighbors[1], &g.nodes[1]);
    ASSERT_EQ(g.nodes[5].neighbors[0], nullptr);
    std::set<size_t> expected_obs = {2};
    ASSERT_EQ(g.negative_weight_observables_set, expected_obs);
    std::set<size_t> expected_detection_events = {0};
    ASSERT_EQ(g.negative_weight_detection_events_set, expected_detection_events);
}

TEST(Graph, AddEdgeWithImpliedWeights) {
    pm::MatchingGraph g(4, 64);
    std::vector<pm::ImpliedWeightUnconverted> implied_weights = {{2, 3, 5}};
    g.add_edge(0, 1, 10, {}, implied_weights);

    ASSERT_EQ(g.edges_to_implied_weights_unconverted.size(), 2);
    ASSERT_TRUE(g.edges_to_implied_weights_unconverted.count(0));
    ASSERT_TRUE(g.edges_to_implied_weights_unconverted.count(1));
    ASSERT_FALSE(g.edges_to_implied_weights_unconverted.count(2));

    ASSERT_EQ(g.edges_to_implied_weights_unconverted[0].size(), 1);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[0][0].size(), 1);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[0][0][0].node1, 2);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[0][0][0].node2, 3);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[0][0][0].implied_weight, 5);

    ASSERT_EQ(g.edges_to_implied_weights_unconverted[1].size(), 1);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[1][0].size(), 1);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[1][0][0].node1, 2);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[1][0][0].node2, 3);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[1][0][0].implied_weight, 5);
}

TEST(Graph, AddBoundaryEdgeWithImpliedWeights) {
    pm::MatchingGraph g(4, 64);
    std::vector<pm::ImpliedWeightUnconverted> implied_weights = {{1, 2, 7}};
    g.add_boundary_edge(0, 10, {}, implied_weights);

    ASSERT_EQ(g.edges_to_implied_weights_unconverted.size(), 1);
    ASSERT_TRUE(g.edges_to_implied_weights_unconverted.count(0));

    ASSERT_EQ(g.edges_to_implied_weights_unconverted[0].size(), 1);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[0][0].size(), 1);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[0][0][0].node1, 1);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[0][0][0].node2, 2);
    ASSERT_EQ(g.edges_to_implied_weights_unconverted[0][0][0].implied_weight, 7);
}

TEST(Graph, ApplyTempReweights) {
    pm::MatchingGraph g(4, 64);
    g.add_edge(0, 1, 10, {0}, {});
    g.add_edge(1, 2, -20, {1}, {});
    
    ASSERT_EQ(g.nodes[0].neighbor_weights[0], 10);
    ASSERT_EQ(g.nodes[1].neighbor_weights[0], 10);
    ASSERT_EQ(g.nodes[1].neighbor_weights[1], 20);
    ASSERT_EQ(g.negative_weight_sum, -20);
    ASSERT_EQ(g.negative_weight_sum_delta, 0);

    // Reweight positive to positive
    std::vector<std::tuple<size_t, int64_t, double>> reweights;
    g.normalising_constant = 1.0;
    // 30.0 * (1.0/2) = 15.0. 15 * 2 = 30.
    reweights.emplace_back(0, 1, 30.0);
    g.apply_temp_reweights(reweights);
    ASSERT_EQ(g.nodes[0].neighbor_weights[0], 30);
    ASSERT_EQ(g.negative_weight_sum, -20); // No change
    ASSERT_EQ(g.negative_weight_sum_delta, 0);

    g.undo_reweights();
    ASSERT_EQ(g.nodes[0].neighbor_weights[0], 10);
    ASSERT_EQ(g.negative_weight_sum_delta, 0);

    // Reweight negative to negative
    // Orig -20. Reweight to -40.0.
    // -40 * 0.5 = -20. -20 * 2 = -40.
    reweights.clear();
    reweights.emplace_back(1, 2, -40.0);
    g.apply_temp_reweights(reweights);
    ASSERT_EQ(g.nodes[1].neighbor_weights[1], 40);
    // Delta: -40 - (-20) = -20.
    // Sum: -20 + (-20) = -40.
    ASSERT_EQ(g.negative_weight_sum, -40);
    ASSERT_EQ(g.negative_weight_sum_delta, -20);

    g.undo_reweights();
    ASSERT_EQ(g.nodes[1].neighbor_weights[1], 20);
    ASSERT_EQ(g.negative_weight_sum, -20);
    ASSERT_EQ(g.negative_weight_sum_delta, 0);

    // Sign flip (should throw)
    reweights.clear();
    reweights.emplace_back(0, 1, -5.0);
    ASSERT_THROW(g.apply_temp_reweights(reweights), std::invalid_argument);
}
