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

#include "pymatching/sparse_blossom/driver/user_graph.h"

#include <cmath>
#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"

double merge_weights_via_probabilities(double a, double b) {
    double p_a = 1 / (1 + std::exp(a));
    double p_b = 1 / (1 + std::exp(b));
    double p_both = p_a * (1 - p_b) + p_b * (1 - p_a);
    return std::log((1 - p_both) / p_both);
}

TEST(StimIO, MergeWeights) {
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(10, 11), pm::merge_weights(10, 11));
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(-4.1, 5), pm::merge_weights(-4.1, 5));
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(25, -24), pm::merge_weights(25, -24));
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(-1, -2), pm::merge_weights(-1, -2));
    ASSERT_FLOAT_EQ(pm::merge_weights(1000, 0), 0);
}

TEST(UserGraph, ConstructGraph) {
    pm::UserGraph graph;
    graph.add_or_merge_boundary_edge(0, {2}, 4.1, 0.1);
    graph.add_or_merge_boundary_edge(0, {3}, 1.0, 0.46, pm::INDEPENDENT);
    graph.add_or_merge_edge(0, 1, {0}, 2.5, 0.4);
    graph.add_or_merge_edge(0, 1, {0}, 2.1, 0.45, pm::INDEPENDENT);
    graph.add_or_merge_edge(1, 2, {1}, -3.5, 0.8);
    graph.add_or_merge_edge(2, 3, {3}, 1.8, 0.3);
    graph.add_or_merge_edge(2, 4, {4}, 2.0, 0.25);
    graph.add_or_merge_boundary_edge(2, {4}, 2.2, 0.25);
    graph.set_boundary({3, 4});
    ASSERT_EQ(graph.get_num_observables(), 5);
    ASSERT_EQ(graph.nodes.size(), 5);
    ASSERT_EQ(graph.nodes[0].neighbors[0].edge_it->node2, SIZE_MAX);
    ASSERT_EQ(graph.nodes[0].neighbors[0].edge_it->weight, pm::merge_weights(4.1, 1.0));
    ASSERT_EQ(graph.nodes[0].neighbors[0].edge_it->error_probability, 0.1 * (1 - 0.46) + 0.46 * (1 - 0.1));
    ASSERT_EQ(graph.nodes[0].neighbors[1].edge_it->node2, 1);
    ASSERT_EQ(graph.nodes[0].neighbors[1].edge_it->weight, pm::merge_weights(2.5, 2.1));
    ASSERT_EQ(graph.nodes[0].neighbors[1].edge_it->error_probability, 0.4 * (1 - 0.45) + 0.45 * (1 - 0.4));
    ASSERT_EQ(graph.nodes[1].neighbors[0].edge_it->weight, pm::merge_weights(2.5, 2.1));
    ASSERT_EQ(graph.nodes[1].neighbors[0].edge_it->error_probability, 0.4 * (1 - 0.45) + 0.45 * (1 - 0.4));
    ASSERT_EQ(graph.nodes[1].neighbors[1].edge_it->weight, -3.5);
    ASSERT_EQ(graph.nodes[1].neighbors[0].edge_it->node1, 0);
    ASSERT_EQ(graph.nodes[2].neighbors[0].edge_it->node1, 1);
    ASSERT_EQ(graph.nodes[2].neighbors[1].edge_it->node2, 3);
    ASSERT_EQ(graph.nodes[0].index_of_neighbor(1), 1);
    ASSERT_EQ(graph.nodes[0].index_of_neighbor(SIZE_MAX), 0);
    ASSERT_EQ(graph.nodes[0].index_of_neighbor(3), SIZE_MAX);
    auto& mwpm = graph.get_mwpm();
    auto& g2 = mwpm.flooder.graph;
    ASSERT_EQ(g2.nodes.size(), 5);
    ASSERT_EQ(g2.nodes[0].neighbors[0], nullptr);
    ASSERT_EQ(
        g2.nodes[0].neighbor_weights[0],
        2 * (pm::weight_int)round(pm::merge_weights(4.1, 1.0) * mwpm.flooder.graph.normalising_constant / 2));
    ASSERT_EQ(g2.nodes[0].neighbor_observables[0], 1 << 2);
    ASSERT_EQ(g2.nodes[0].neighbors[1], &g2.nodes[1]);
    ASSERT_EQ(
        g2.nodes[0].neighbor_weights[1],
        2 * (pm::weight_int)round(pm::merge_weights(2.5, 2.1) * mwpm.flooder.graph.normalising_constant / 2));
    ASSERT_EQ(g2.nodes[0].neighbor_observables[1], 1 << 0);
    ASSERT_EQ(
        g2.nodes[1].neighbor_weights[1], 2 * (pm::weight_int)round(3.5 * mwpm.flooder.graph.normalising_constant / 2));
    ASSERT_EQ(g2.nodes[2].neighbors[0], nullptr);
    ASSERT_EQ(
        g2.nodes[2].neighbor_weights[0], 2 * (pm::weight_int)round(1.8 * mwpm.flooder.graph.normalising_constant / 2));
    ASSERT_EQ(g2.nodes[2].neighbor_observables[0], 1 << 3);
    ASSERT_EQ(g2.nodes[2].neighbors[1], &g2.nodes[1]);
    ASSERT_EQ(
        g2.nodes[2].neighbor_weights[1], 2 * (pm::weight_int)round(3.5 * mwpm.flooder.graph.normalising_constant / 2));
    ASSERT_EQ(g2.nodes[2].neighbor_observables[1], 1 << 1);
    ASSERT_EQ(g2.nodes[2].neighbors.size(), 2);
    ASSERT_EQ(g2.nodes[3].neighbors.size(), 0);
    ASSERT_EQ(
        mwpm.flooder.negative_weight_sum,
        -2 * (pm::total_weight_int)round(3.5 * mwpm.flooder.graph.normalising_constant / 2));
    std::set<size_t> dets_exp = {1, 2};
    ASSERT_EQ(mwpm.flooder.graph.negative_weight_detection_events_set.size(), 2);
    std::set<size_t> obs_exp = {1};
    ASSERT_EQ(mwpm.flooder.graph.negative_weight_observables_set, obs_exp);
    ASSERT_EQ(mwpm.flooder.negative_weight_obs_mask, 1 << 1);
    std::vector<uint64_t> dets_exp_vec = {1, 2};
    ASSERT_EQ(mwpm.flooder.negative_weight_detection_events, dets_exp_vec);
    std::vector<size_t> obs_exp_vec = {1};
    ASSERT_EQ(mwpm.flooder.negative_weight_observables, obs_exp_vec);
}

TEST(UserGraph, AddNoise) {
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

TEST(UserGraph, NodesAlongShortestPath) {
    pm::UserGraph graph;
    graph.add_or_merge_boundary_edge(0, {0}, 1, -1);
    graph.add_or_merge_edge(0, 1, {1}, 1, -1);
    graph.add_or_merge_edge(1, 2, {2}, 1, -1);
    graph.add_or_merge_edge(2, 3, {3}, 1, -1);
    graph.add_or_merge_edge(3, 4, {4}, 1, -1);
    graph.add_or_merge_edge(4, 5, {5}, 1, -1);
    graph.set_boundary({5});

    {
        std::vector<size_t> nodes;
        graph.get_nodes_on_shortest_path_from_source(4, 0, nodes);
        std::vector<size_t> nodes_expected = {4, 3, 2, 1, 0};
        ASSERT_EQ(nodes, nodes_expected);
    }

    {
        std::vector<size_t> nodes;
        graph.get_nodes_on_shortest_path_from_source(1, SIZE_MAX, nodes);
        std::vector<size_t> nodes_expected = {1, 0};
        ASSERT_EQ(nodes, nodes_expected);
    }

    {
        std::vector<size_t> nodes;
        graph.get_nodes_on_shortest_path_from_source(SIZE_MAX, 3, nodes);
        std::vector<size_t> nodes_expected = {4, 3};
        ASSERT_EQ(nodes, nodes_expected);
    }

    {
        std::vector<size_t> nodes;
        graph.get_nodes_on_shortest_path_from_source(5, 3, nodes);
        std::vector<size_t> nodes_expected = {4, 3};
        ASSERT_EQ(nodes, nodes_expected);
    }
}

TEST(UserGraph, DecodeUserGraphDetectionEventOnBoundaryNode) {
    {
        pm::UserGraph graph;
        graph.add_or_merge_edge(0, 1, {0}, 1.0, -1);
        graph.add_or_merge_edge(1, 2, {1}, 1.0, -1);
        graph.set_boundary({2});
        auto& mwpm = graph.get_mwpm();
        pm::ExtendedMatchingResult res(mwpm.flooder.graph.num_observables);
        pm::decode_detection_events(mwpm, {2}, res.obs_crossed.data(), res.weight);
    }

    {
        pm::UserGraph graph;
        graph.add_or_merge_edge(0, 1, {0}, -1.0, -1);
        graph.add_or_merge_edge(1, 2, {1}, 1.0, -1);
        graph.set_boundary({2});
        auto& mwpm = graph.get_mwpm();
        pm::ExtendedMatchingResult res(mwpm.flooder.graph.num_observables);
        pm::decode_detection_events(mwpm, {2}, res.obs_crossed.data(), res.weight);
    }
}

/// Tests for iter_dem_instructions_include_correlations

/// A helper struct to capture the data passed to the handler function for later validation in tests.
struct HandledError {
    double probability;
    size_t node1;
    size_t node2;
    std::vector<size_t> observables;

    // Equality operator for easy comparison with ASSERT_EQ.
    bool operator==(const HandledError& other) const {
        return probability == other.probability && node1 == other.node1 && node2 == other.node2 &&
               observables == other.observables;
    }
};

struct TestHandler {
    mutable std::vector<HandledError> handled_errors;

    void operator()(double p, const std::vector<size_t> dets, const std::vector<size_t>& observables) const {
        if (dets.size() == 1) {
            handled_errors.push_back({p, dets[0], SIZE_MAX, observables});
        } else if (dets.size() == 2) {
            handled_errors.push_back({p, dets[0], dets[1], observables});
        }
    }
};

// Test with an empty detector error model. The handler should not be called.
TEST(IterDemInstructionsTest, EmptyDem) {
    stim::DetectorErrorModel dem;
    TestHandler handler;
    pm::iter_dem_instructions_include_correlations(dem, handler);
    ASSERT_TRUE(handler.handled_errors.empty());
}

// Test a simple error involving one detector, which is an error on the boundary.
TEST(IterDemInstructionsTest, SingleDetectorErrorToBoundary) {
    stim::DetectorErrorModel dem("error(0.1) D0");
    TestHandler handler;
    pm::iter_dem_instructions_include_correlations(dem, handler);
    ASSERT_EQ(handler.handled_errors.size(), 1);
    EXPECT_EQ(handler.handled_errors[0], (HandledError{0.1, 0, SIZE_MAX, {}}));
}

// Test a standard error between two detectors.
TEST(IterDemInstructionsTest, TwoDetectorError) {
    stim::DetectorErrorModel dem("error(0.25) D5 D10");
    TestHandler handler;
    pm::iter_dem_instructions_include_correlations(dem, handler);
    ASSERT_EQ(handler.handled_errors.size(), 1);
    EXPECT_EQ(handler.handled_errors[0], (HandledError{0.25, 5, 10, {}}));
}

// Test an error that also flips a logical observable.
TEST(IterDemInstructionsTest, ErrorWithOneObservable) {
    stim::DetectorErrorModel dem("error(0.125) D1 D2 L0");
    TestHandler handler;
    pm::iter_dem_instructions_include_correlations(dem, handler);
    ASSERT_EQ(handler.handled_errors.size(), 1);
    EXPECT_EQ(handler.handled_errors[0], (HandledError{0.125, 1, 2, {0}}));
}

// Test an error that flips multiple logical observables.
TEST(IterDemInstructionsTest, ErrorWithMultipleObservables) {
    stim::DetectorErrorModel dem("error(0.3) D3 D4 L1 L3");
    TestHandler handler;
    pm::iter_dem_instructions_include_correlations(dem, handler);
    ASSERT_EQ(handler.handled_errors.size(), 1);
    EXPECT_EQ(handler.handled_errors[0], (HandledError{0.3, 3, 4, {1, 3}}));
}

// Test an error with probability 0. It should be ignored.
TEST(IterDemInstructionsTest, ZeroProbabilityError) {
    stim::DetectorErrorModel dem("error(0.0) D0 D1");
    TestHandler handler;
    pm::iter_dem_instructions_include_correlations(dem, handler);
    ASSERT_TRUE(handler.handled_errors.empty());
}

// Test an error involving more than two detectors (a hyperedge). This should be ignored.
TEST(IterDemInstructionsTest, ThreeDetectorErrorIsIgnored) {
    stim::DetectorErrorModel dem("error(0.1) D0 D1 D2");
    TestHandler handler;
    pm::iter_dem_instructions_include_correlations(dem, handler);
    ASSERT_TRUE(handler.handled_errors.empty());
}

// Test a decomposed error instruction. The handler should be called for each component.
TEST(IterDemInstructionsTest, DecomposedError) {
    stim::DetectorErrorModel dem("error(0.1) D0 D1 ^ D2 D3 L0 ^ D4");
    TestHandler handler;
    pm::iter_dem_instructions_include_correlations(dem, handler);
    ASSERT_EQ(handler.handled_errors.size(), 3);

    // std::vector<HandledError> expected = {{0.1, 0, 1, {}}, {0.1, 2, 3, {0}}, {0.1, 4, SIZE_MAX, {}}};
    EXPECT_EQ(handler.handled_errors[0], (HandledError{0.1, 0, 1, {}}));
    EXPECT_EQ(handler.handled_errors[1], (HandledError{0.1, 2, 3, {0}}));
    EXPECT_EQ(handler.handled_errors[2], (HandledError{0.1, 4, SIZE_MAX, {}}));
}

// Test a decomposed error where one of the components is a hyperedge and should be ignored.
TEST(IterDemInstructionsTest, DecomposedErrorWithIgnoredComponent) {
    stim::DetectorErrorModel dem("error(0.15) D0 D1 ^ D2 D3 D4 ^ D5 D6 L2");
    TestHandler handler;
    pm::iter_dem_instructions_include_correlations(dem, handler);
    ASSERT_EQ(handler.handled_errors.size(), 2);

    std::vector<HandledError> expected = {{0.15, 0, 1, {}}, {0.15, 5, 6, {2}}};
    EXPECT_EQ(handler.handled_errors, expected);
}

// Test a complex DEM with multiple instruction types and edge cases combined.
TEST(IterDemInstructionsTest, CombinedComplexDem) {
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0            # Instruction 1: Simple
        error(0.2) D1 D2 L0      # Instruction 2: Two detectors, one observable
        error(0.3) D3 D4 D5 ^ D6 # Instruction 3: Hyperedge ignored, second component handled
        error(0.0) D7            # Instruction 4: Zero probability, ignored
        error(0.4) D8 ^ D9 L1    # Instruction 5: Decomposed
    )DEM");
    TestHandler handler;
    pm::iter_dem_instructions_include_correlations(dem, handler);

    // We expect 5 errors to be handled from this model.
    // 1. {0.1, 0, MAX, {}} from "error(0.1) D0"
    // 2. {0.2, 1, 2, {0}} from "error(0.2) D1 D2 L0"
    // 3. {0.3, 6, MAX, {}} from the second part of "error(0.3) D3 D4 D5 ^ D6" (first part is ignored)
    // 4. {0.4, 8, MAX, {}} from the first part of "error(0.4) D8 ^ D9 L1"
    // 5. {0.4, 9, MAX, {1}} from the second part of "error(0.4) D8 ^ D9 L1"
    ASSERT_EQ(handler.handled_errors.size(), 5);

    std::vector<HandledError> expected = {
        {0.1, 0, SIZE_MAX, {}},
        {0.2, 1, 2, {0}},
        {0.3, 6, SIZE_MAX, {}},
        {0.4, 8, SIZE_MAX, {}},
        {0.4, 9, SIZE_MAX, {1}}};

    // stim's iter_flatten_error_instructions processes instructions sequentially,
    // so the order of handled errors should be deterministic.
    EXPECT_EQ(handler.handled_errors, expected);
}
