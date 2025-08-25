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
#include <limits>
#include <stdexcept>

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
        pm::decode_detection_events(mwpm, {2}, res.obs_crossed.data(), res.weight, /*enable_correlations=*/false);
    }

    {
        pm::UserGraph graph;
        graph.add_or_merge_edge(0, 1, {0}, -1.0, -1);
        graph.add_or_merge_edge(1, 2, {1}, 1.0, -1);
        graph.set_boundary({2});
        auto& mwpm = graph.get_mwpm();
        pm::ExtendedMatchingResult res(mwpm.flooder.graph.num_observables);
        pm::decode_detection_events(mwpm, {2}, res.obs_crossed.data(), res.weight, /*enable_correlations=*/false);
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
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);
    ASSERT_TRUE(handler.handled_errors.empty());
    ASSERT_TRUE(joint_probabilities.empty());
}

// Test a simple error involving one detector, which is an error on the boundary.
TEST(IterDemInstructionsTest, SingleDetectorErrorToBoundary) {
    stim::DetectorErrorModel dem("error(0.1) D0");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);

    // Check handler calls
    ASSERT_EQ(handler.handled_errors.size(), 1);
    EXPECT_EQ(handler.handled_errors[0], (HandledError{0.1, 0, SIZE_MAX, {}}));

    // Check joint probabilities (marginal probability in this case)
    std::pair<size_t, size_t> key = {0, SIZE_MAX};
    ASSERT_EQ(joint_probabilities.size(), 1);
    ASSERT_EQ(joint_probabilities[key].size(), 1);
    EXPECT_DOUBLE_EQ(joint_probabilities[key][key], 0.1);
}

// Test a standard error between two detectors.
TEST(IterDemInstructionsTest, TwoDetectorError) {
    stim::DetectorErrorModel dem("error(0.25) D5 D10");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);

    ASSERT_EQ(handler.handled_errors.size(), 1);
    EXPECT_EQ(handler.handled_errors[0], (HandledError{0.25, 5, 10, {}}));

    std::pair<size_t, size_t> key = {5, 10};
    EXPECT_DOUBLE_EQ(joint_probabilities[key][key], 0.25);
}

// Test a standard error between two detectors where they are not sorted in the DEM.
TEST(IterDemInstructionsTest, TwoDetectorErrorNotSorted) {
    stim::DetectorErrorModel dem("error(0.25) D10 D5");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);

    ASSERT_EQ(handler.handled_errors.size(), 1);
    EXPECT_EQ(handler.handled_errors[0], (HandledError{0.25, 5, 10, {}}));

    std::pair<size_t, size_t> key = {5, 10};
    EXPECT_DOUBLE_EQ(joint_probabilities[key][key], 0.25);
}

// Test an error that also flips a logical observable.
TEST(IterDemInstructionsTest, ErrorWithOneObservable) {
    stim::DetectorErrorModel dem("error(0.125) D1 D2 L0");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);

    ASSERT_EQ(handler.handled_errors.size(), 1);
    EXPECT_EQ(handler.handled_errors[0], (HandledError{0.125, 1, 2, {0}}));

    std::pair<size_t, size_t> key = {1, 2};
    EXPECT_DOUBLE_EQ(joint_probabilities[key][key], 0.125);
}

TEST(IterDemInstructionsTest, ErrorWithMultipleObservables) {
    stim::DetectorErrorModel dem("error(0.3) D3 D4 L1 L3");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);

    ASSERT_EQ(handler.handled_errors.size(), 1);
    EXPECT_EQ(handler.handled_errors[0], (HandledError{0.3, 3, 4, {1, 3}}));

    std::pair<size_t, size_t> key = {3, 4};
    EXPECT_DOUBLE_EQ(joint_probabilities[key][key], 0.3);
}

TEST(IterDemInstructionsTest, ZeroProbabilityError) {
    stim::DetectorErrorModel dem("error(0.0) D0 D1");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);
    ASSERT_TRUE(handler.handled_errors.empty());
    ASSERT_TRUE(joint_probabilities.empty());
}

TEST(IterDemInstructionsTest, ThreeDetectorErrorThrowsInvalidArgument) {
    stim::DetectorErrorModel dem("error(0.1) D0 D1 D2");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    ASSERT_THROW(
        pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities), std::invalid_argument);
}

// Test a decomposed error instruction. The handler should be called for each component.
TEST(IterDemInstructionsTest, DecomposedError) {
    stim::DetectorErrorModel dem("error(0.1) D0 D1 ^ D2 D3 L0 ^ D4");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);

    // Check handler calls
    ASSERT_EQ(handler.handled_errors.size(), 3);
    std::vector<HandledError> expected_handled = {{0.1, 0, 1, {}}, {0.1, 2, 3, {0}}, {0.1, 4, SIZE_MAX, {}}};
    EXPECT_EQ(handler.handled_errors, expected_handled);

    // Check joint probabilities
    std::pair<size_t, size_t> key01 = {0, 1};
    std::pair<size_t, size_t> key23 = {2, 3};
    std::pair<size_t, size_t> key4B = {4, SIZE_MAX};

    // Marginal probabilities
    EXPECT_DOUBLE_EQ(joint_probabilities[key01][key01], 0.1);
    EXPECT_DOUBLE_EQ(joint_probabilities[key23][key23], 0.1);
    EXPECT_DOUBLE_EQ(joint_probabilities[key4B][key4B], 0.1);

    // Joint probabilities between components
    EXPECT_DOUBLE_EQ(joint_probabilities[key01][key23], 0.1);
    EXPECT_DOUBLE_EQ(joint_probabilities[key23][key01], 0.1);  // Symmetric
    EXPECT_DOUBLE_EQ(joint_probabilities[key01][key4B], 0.1);
    EXPECT_DOUBLE_EQ(joint_probabilities[key4B][key01], 0.1);  // Symmetric
    EXPECT_DOUBLE_EQ(joint_probabilities[key23][key4B], 0.1);
    EXPECT_DOUBLE_EQ(joint_probabilities[key4B][key23], 0.1);  // Symmetric
}

// Test that a decomposed error with a hyperedge component throws an exception.
TEST(IterDemInstructionsTest, DecomposedErrorWithHyperedgeThrows) {
    stim::DetectorErrorModel dem("error(0.15) D0 D1 ^ D2 D3 D4 ^ D5 D6 L2");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;

    // Assert that the function throws std::invalid_argument when processing the DEM.
    ASSERT_THROW(
        pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities), std::invalid_argument);
}

// Test that a decomposed error with an undetectable component throws an exception.
TEST(IterDemInstructionsTest, DecomposedErrorWithUndetectableErrorThrows) {
    stim::DetectorErrorModel dem("error(0.15) L0 ^ D2 D4 ^ D5 D6 L2");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;

    // Assert that the function throws std::invalid_argument when processing the DEM.
    ASSERT_THROW(
        pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities), std::invalid_argument);

    stim::DetectorErrorModel dem2("error(0.15) D2 D4 ^ D5 D6 L2 ^ L1");
    // Assert that the function throws std::invalid_argument when processing the DEM.
    ASSERT_THROW(
        pm::iter_dem_instructions_include_correlations(dem2, handler, joint_probabilities), std::invalid_argument);
}

// Test a complex DEM with multiple instruction types and edge cases combined.
TEST(IterDemInstructionsTest, CombinedComplexDem) {
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0            # Instruction 1: Simple
        error(0.3) L0            # Instruction 2: Undetectable error, ignored
        error(0.2) D1 D2 L0      # Instruction 2: Two detectors, one observable
        error(0.0) D7            # Instruction 3: Zero probability, ignored
        error(0.4) D8 ^ D9 L1    # Instruction 4: Decomposed
    )DEM");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);

    ASSERT_EQ(handler.handled_errors.size(), 4);
    
    std::vector<HandledError> expected = {
        {0.1, 0, SIZE_MAX, {}}, {0.2, 1, 2, {0}}, {0.4, 8, SIZE_MAX, {}}, {0.4, 9, SIZE_MAX, {1}}};
    EXPECT_EQ(handler.handled_errors, expected);

    // Check joint probabilities
    std::pair<size_t, size_t> key0B = {0, SIZE_MAX};
    std::pair<size_t, size_t> key12 = {1, 2};
    std::pair<size_t, size_t> key8B = {8, SIZE_MAX};
    std::pair<size_t, size_t> key9B = {9, SIZE_MAX};

    // Marginal probabilities from each instruction
    EXPECT_DOUBLE_EQ(joint_probabilities[key0B][key0B], 0.1);
    EXPECT_DOUBLE_EQ(joint_probabilities[key12][key12], 0.2);
    EXPECT_DOUBLE_EQ(joint_probabilities[key8B][key8B], 0.4);
    EXPECT_DOUBLE_EQ(joint_probabilities[key9B][key9B], 0.4);

    // Joint probability from the last instruction
    EXPECT_DOUBLE_EQ(joint_probabilities[key8B][key9B], 0.4);
    EXPECT_DOUBLE_EQ(joint_probabilities[key9B][key8B], 0.4);

    // Check that there are no other joint probabilities
    EXPECT_EQ(joint_probabilities[key0B].count(key12), 0);
}

double bernoulli_xor(double p1, double p2) {
    return p1 * (1 - p2) + p2 * (1 - p1);
}

TEST(IterDemInstructionsTest, MoreThanEightComponents) {
    stim::DetectorErrorModel dem("error(0.1) D0 ^ D1 ^ D2 ^ D3 ^ D4 ^ D5 ^ D6 ^ D7 ^ D8");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);
}

// Tests that multiple error instructions on the same edge correctly combine their probabilities.
TEST(IterDemInstructionsTest, MultipleErrorsOnSameEdgeCombine) {
    stim::DetectorErrorModel dem("error(0.1) D0 D1\n error(0.2) D0 D1");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);

    // Expected probability = 0.1*(1-0.2) + 0.2*(1-0.1) = 0.08 + 0.18 = 0.26
    double expected_p = bernoulli_xor(0.1, 0.2);

    std::pair<size_t, size_t> key = {0, 1};
    EXPECT_DOUBLE_EQ(joint_probabilities[key][key], expected_p);
}

// Tests how marginal and joint probabilities are combined across different decomposed error instructions.
TEST(IterDemInstructionsTest, ComplexCombinationOfErrors) {
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0 ^ D1
        error(0.2) D0 ^ D2
    )DEM");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities);

    std::pair<size_t, size_t> k0 = {0, SIZE_MAX};
    std::pair<size_t, size_t> k1 = {1, SIZE_MAX};
    std::pair<size_t, size_t> k2 = {2, SIZE_MAX};

    // Marginal probabilities
    EXPECT_DOUBLE_EQ(joint_probabilities[k0][k0], bernoulli_xor(0.1, 0.2));  // 0.26
    EXPECT_DOUBLE_EQ(joint_probabilities[k1][k1], 0.1);
    EXPECT_DOUBLE_EQ(joint_probabilities[k2][k2], 0.2);

    // Joint probabilities
    EXPECT_DOUBLE_EQ(joint_probabilities[k0][k1], 0.1);
    EXPECT_DOUBLE_EQ(joint_probabilities[k1][k0], 0.1);
    EXPECT_DOUBLE_EQ(joint_probabilities[k0][k2], 0.2);
    EXPECT_DOUBLE_EQ(joint_probabilities[k2][k0], 0.2);

    // No instruction connects D1 and D2, so their joint probability should be 0.
    EXPECT_DOUBLE_EQ(joint_probabilities[k1][k2], 0.0);
}

// Test that an error greater than 0.5 results in a throw.
TEST(IterDemInstructionsTest, ProbabilityGreaterThanHalfThrows) {
    stim::DetectorErrorModel dem("error(0.51) D0 D2");
    TestHandler handler;
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    ASSERT_THROW(
        pm::iter_dem_instructions_include_correlations(dem, handler, joint_probabilities), std::invalid_argument);
}

TEST(UserGraph, PopulateImpliedEdgeWeights) {
    pm::UserGraph graph;
    graph.add_or_merge_edge(0, 1, {}, 0.0, 0.26);
    graph.add_or_merge_edge(2, 3, {}, 0.0, 0.1);

    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>> joint_probabilities;
    joint_probabilities[{0, 1}][{0, 1}] = 0.26;
    joint_probabilities[{0, 1}][{2, 3}] = 0.1;
    joint_probabilities[{2, 3}][{0, 1}] = 0.1;
    joint_probabilities[{2, 3}][{2, 3}] = 0.1;

    graph.populate_implied_edge_weights(joint_probabilities);

    auto it_01 = std::find_if(graph.edges.begin(), graph.edges.end(), [](const pm::UserEdge& edge) {
        return edge.node1 == 0 && edge.node2 == 1;
    });
    ASSERT_NE(it_01, graph.edges.end());
    ASSERT_EQ(it_01->implied_weights_for_other_edges.size(), 1);
    const auto& implied_01 = it_01->implied_weights_for_other_edges[0];
    ASSERT_EQ(implied_01.node1, 2);
    ASSERT_EQ(implied_01.node2, 3);

    double p_01 = 0.1 / 0.26;
    double w_01 = pm::to_weight_for_correlations(p_01);
    ASSERT_EQ(implied_01.implied_weight, w_01);

    auto it_23 = std::find_if(graph.edges.begin(), graph.edges.end(), [](const pm::UserEdge& edge) {
        return edge.node1 == 2 && edge.node2 == 3;
    });

    ASSERT_NE(it_23, graph.edges.end());
    ASSERT_EQ(it_23->implied_weights_for_other_edges.size(), 1);
    const auto& implied_23 = it_23->implied_weights_for_other_edges[0];
    ASSERT_EQ(implied_23.node1, 0);
    ASSERT_EQ(implied_23.node2, 1);
    ASSERT_NEAR(implied_23.implied_weight, 0.0, 0.00001);
}

TEST(UserGraph, ConvertImpliedWeights) {
    pm::UserGraph user_graph;
    user_graph.add_or_merge_edge(0, 1, {}, 1.0, 0.1);
    user_graph.add_or_merge_edge(2, 3, {}, 1.0, 0.1);
    user_graph.add_or_merge_boundary_edge(4, {}, 1.0, 0.1);

    auto& edge1 = *std::find_if(user_graph.edges.begin(), user_graph.edges.end(), [](const pm::UserEdge& e) {
        return e.node1 == 0;
    });
    edge1.implied_weights_for_other_edges.push_back({2, 3, 5});
    edge1.implied_weights_for_other_edges.push_back({4, SIZE_MAX, 7});

    pm::MatchingGraph matching_graph = user_graph.to_matching_graph(100);

    auto& node0_neighbors = matching_graph.nodes[0].neighbors;
    auto it0 = std::find(node0_neighbors.begin(), node0_neighbors.end(), &matching_graph.nodes[1]);
    size_t index_of_1_in_0 = std::distance(node0_neighbors.begin(), it0);

    auto& node1_neighbors = matching_graph.nodes[1].neighbors;
    auto it1 = std::find(node1_neighbors.begin(), node1_neighbors.end(), &matching_graph.nodes[0]);
    size_t index_of_0_in_1 = std::distance(node1_neighbors.begin(), it1);

    auto& node2_neighbors = matching_graph.nodes[2].neighbors;
    auto it2 = std::find(node2_neighbors.begin(), node2_neighbors.end(), &matching_graph.nodes[3]);
    size_t index_of_3_in_2 = std::distance(node2_neighbors.begin(), it2);

    auto& node3_neighbors = matching_graph.nodes[3].neighbors;
    auto it3 = std::find(node3_neighbors.begin(), node3_neighbors.end(), &matching_graph.nodes[2]);
    size_t index_of_2_in_3 = std::distance(node3_neighbors.begin(), it3);

    auto& node4_neighbors = matching_graph.nodes[4].neighbors;
    auto it4 = std::find(node4_neighbors.begin(), node4_neighbors.end(), nullptr);
    size_t index_of_boundary_in_4 = std::distance(node4_neighbors.begin(), it4);

    auto& implied_weights = matching_graph.nodes[0].neighbor_implied_weights[index_of_1_in_0];
    ASSERT_EQ(implied_weights.size(), 2);
    ASSERT_EQ(implied_weights[0].edge0_ptr, &matching_graph.nodes[2].neighbor_weights[index_of_3_in_2]);
    ASSERT_EQ(implied_weights[0].edge1_ptr, &matching_graph.nodes[3].neighbor_weights[index_of_2_in_3]);
    ASSERT_EQ(implied_weights[0].implied_weight, 10);
    ASSERT_EQ(implied_weights[1].edge0_ptr, &matching_graph.nodes[4].neighbor_weights[index_of_boundary_in_4]);
    ASSERT_EQ(implied_weights[1].edge1_ptr, nullptr);
    ASSERT_EQ(implied_weights[1].implied_weight, 14);

    auto& implied_weights_rev = matching_graph.nodes[1].neighbor_implied_weights[index_of_0_in_1];
    ASSERT_EQ(implied_weights_rev.size(), 2);
}

TEST(UserGraph, ConvertImpliedWeights_NoRules) {
    pm::UserGraph user_graph;
    user_graph.add_or_merge_edge(0, 1, {}, 1.0, 0.1);
    user_graph.add_or_merge_edge(2, 3, {}, 1.0, 0.1);
    user_graph.add_or_merge_boundary_edge(4, {}, 1.0, 0.1);

    pm::MatchingGraph matching_graph = user_graph.to_matching_graph(100);

    for (const auto& node : matching_graph.nodes) {
        for (const auto& implied_weights_vec : node.neighbor_implied_weights) {
            ASSERT_TRUE(implied_weights_vec.empty());
        }
    }
}

TEST(UserGraph, ConvertImpliedWeights_EmptyRules) {
    pm::UserGraph user_graph;
    user_graph.add_or_merge_edge(0, 1, {}, 1.0, 0.1);
    user_graph.add_or_merge_edge(2, 3, {}, 1.0, 0.1);

    auto& edge = *std::find_if(user_graph.edges.begin(), user_graph.edges.end(), [](const pm::UserEdge& e) {
        return e.node1 == 0;
    });
    edge.implied_weights_for_other_edges = {};

    pm::MatchingGraph matching_graph = user_graph.to_matching_graph(100);

    for (const auto& node : matching_graph.nodes) {
        for (const auto& implied_weights_vec : node.neighbor_implied_weights) {
            ASSERT_TRUE(implied_weights_vec.empty());
        }
    }
}

TEST(UserGraph, GetEdgeOrBoundaryEdgeWeight) {
    pm::UserGraph user_graph;
    user_graph.add_or_merge_edge(2, 3, {}, 3.5, 0.1);
    user_graph.add_or_merge_boundary_edge(5, {}, 5.1, 0.2);

    double edge_weight = 0;
    bool has_edge = user_graph.get_edge_or_boundary_edge_weight(2, 3, edge_weight);
    ASSERT_EQ(edge_weight, 3.5);
    ASSERT_TRUE(has_edge);
    double boundary_edge_weight = 0;
    bool has_boundary_edge = user_graph.get_edge_or_boundary_edge_weight(5, SIZE_MAX, boundary_edge_weight);
    ASSERT_EQ(boundary_edge_weight, 5.1);
    ASSERT_TRUE(has_boundary_edge);
    double w;
    bool has_non_existent_edge = user_graph.get_edge_or_boundary_edge_weight(2, 4, w);
    ASSERT_FALSE(has_non_existent_edge);
    bool has_non_existent_edge_2 = user_graph.get_edge_or_boundary_edge_weight(10, 0, w);
    ASSERT_FALSE(has_non_existent_edge_2);
}