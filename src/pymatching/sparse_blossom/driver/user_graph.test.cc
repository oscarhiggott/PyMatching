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
#include <cstdint>
#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"

double convert_p_to_double_weight(double p) {
    return std::log((1 - p) / p);
};

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
    ASSERT_EQ(graph.get_num_nodes(), 5);

    // Graph should contain a boundary edge with node 0.
    ASSERT_TRUE(graph.contains(0, SIZE_MAX));

    // Assertions for boundary edge {0, SIZE_MAX}
    std::optional<pm::DetectorEdgeData> edge_data_0_boundary = graph.get_edge_data(0, SIZE_MAX);
    ASSERT_TRUE(edge_data_0_boundary.has_value());
    ASSERT_EQ(edge_data_0_boundary->weight, pm::merge_weights(4.1, 1.0));
    ASSERT_EQ(edge_data_0_boundary->error_probability, 0.1 * (1 - 0.46) + 0.46 * (1 - 0.1));

    // Assertions for {0,1}
    ASSERT_TRUE(graph.contains(0, 1));
    std::optional<pm::DetectorEdgeData> edge_data_0_1 = graph.get_edge_data(0, 1);
    ASSERT_TRUE(edge_data_0_1.has_value());
    ASSERT_EQ(edge_data_0_1->weight, pm::merge_weights(2.5, 2.1));
    ASSERT_EQ(edge_data_0_1->error_probability, 0.4 * (1 - 0.45) + 0.45 * (1 - 0.4));
    //
    //  Assertions for {1,0}
    ASSERT_TRUE(graph.contains(0, 1));
    std::optional<pm::DetectorEdgeData> edge_data_1_0 = graph.get_edge_data(1, 0);
    ASSERT_TRUE(edge_data_1_0.has_value());
    ASSERT_EQ(edge_data_1_0->weight, pm::merge_weights(2.5, 2.1));
    ASSERT_EQ(edge_data_1_0->error_probability, 0.4 * (1 - 0.45) + 0.45 * (1 - 0.4));

    //  Assertions for {1,2}
    ASSERT_TRUE(graph.contains(1, 2));

    std::optional<pm::DetectorEdgeData> edge_data_1_2 = graph.get_edge_data(1, 2);
    ASSERT_TRUE(edge_data_1_2.has_value());
    ASSERT_EQ(edge_data_1_2->weight, -3.5);

    ASSERT_TRUE(graph.contains(2, 1));
    ASSERT_TRUE(graph.contains(2, 3));

    // Assert that a non-inserted edge isn't marked as contained.
    ASSERT_FALSE(graph.contains(100, 100));

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

TEST(ConvertProbabilityToWeightSingleTest, AllCases) {
    EXPECT_EQ(0u, pm::convert_probability_to_weight(0.5));
    EXPECT_EQ(1u, pm::convert_probability_to_weight(0.25));
    EXPECT_EQ(UINT32_MAX, pm::convert_probability_to_weight(0.75));
    EXPECT_EQ(11u, pm::convert_probability_to_weight(0.00001));
    EXPECT_EQ(UINT32_MAX - 10, pm::convert_probability_to_weight(0.99999));
    EXPECT_EQ(0u, pm::convert_probability_to_weight(0.0));
    EXPECT_EQ(0u, pm::convert_probability_to_weight(1.0));
    EXPECT_EQ(0u, pm::convert_probability_to_weight(2.0));
    EXPECT_EQ(0u, pm::convert_probability_to_weight(-1.0));
    EXPECT_EQ(0u, pm::convert_probability_to_weight(std::numeric_limits<double>::quiet_NaN()));
}

TEST(UserGraphFromDem, EmptyDEM) {
    stim::DetectorErrorModel dem;  // Empty model
    pm::UserGraph graph = pm::detector_error_model_to_user_graph(dem);
    ASSERT_EQ(graph.get_num_nodes(), 0);
    ASSERT_EQ(graph.get_num_edges(), 0);
}

TEST(DetectorErrorModelToUserGraph, SingleDetectorError) {
    stim::DetectorErrorModel dem("error(0.1) D0");
    pm::UserGraph graph = pm::detector_error_model_to_user_graph(dem);

    pm::DetectorEdgeId d0(0);

    ASSERT_EQ(graph.edge_map.size(), 1);
    ASSERT_TRUE(graph.edge_map.count(d0));

    const auto& data = graph.edge_map.at(d0);
    EXPECT_EQ(data.detectors, d0);
    EXPECT_TRUE(data.observables.empty());
    EXPECT_DOUBLE_EQ(data.error_probability, 0.1);  // bernoulli_xor(0, 0.1)
    EXPECT_DOUBLE_EQ(data.weight, convert_p_to_double_weight(0.1));

    // From populate_implied_edge_weights
    // conditional_groups[d0][d0] = 0.1
    EXPECT_DOUBLE_EQ(data.correlated_probabilities_sum, 0.1);
    EXPECT_TRUE(data.implied_weights_for_other_edges.empty());
}

TEST(DetectorErrorModelToUserGraph, TwoDetectorError) {
    stim::DetectorErrorModel dem("error(0.2) D1 D3");
    pm::UserGraph graph = pm::detector_error_model_to_user_graph(dem);

    pm::DetectorEdgeId d13(1, 3);

    ASSERT_EQ(graph.edge_map.size(), 1);
    ASSERT_TRUE(graph.edge_map.count(d13));

    const auto& data = graph.edge_map.at(d13);
    EXPECT_EQ(data.detectors, d13);
    EXPECT_TRUE(data.observables.empty());
    EXPECT_DOUBLE_EQ(data.error_probability, 0.2);
    EXPECT_DOUBLE_EQ(data.weight, convert_p_to_double_weight(0.2));

    // conditional_groups[d13][d13] = 0.2
    EXPECT_DOUBLE_EQ(data.correlated_probabilities_sum, 0.2);
    EXPECT_TRUE(data.implied_weights_for_other_edges.empty());
}

TEST(DetectorErrorModelToUserGraph, ThreeDetectorErrorIsIgnored) {
    stim::DetectorErrorModel dem("error(0.3) D0 D1 D2");
    pm::UserGraph graph = pm::detector_error_model_to_user_graph(dem);

    EXPECT_TRUE(graph.edge_map.empty());
}

TEST(DetectorErrorModelToUserGraph, CorrelatedErrorSimple) {
    stim::DetectorErrorModel dem("error(0.6) D0 ^ D1 L0");
    pm::UserGraph graph = pm::detector_error_model_to_user_graph(dem);

    pm::DetectorEdgeId d0(0);
    pm::DetectorEdgeId d1(1);

    ASSERT_EQ(graph.edge_map.size(), 2);  // Edges D0 and D1

    // Check D0
    ASSERT_TRUE(graph.edge_map.count(d0));
    const auto& data_d0 = graph.edge_map.at(d0);
    EXPECT_EQ(data_d0.detectors, d0);
    EXPECT_TRUE(data_d0.observables.empty());
    EXPECT_DOUBLE_EQ(data_d0.error_probability, 0.6);
    EXPECT_DOUBLE_EQ(data_d0.weight, convert_p_to_double_weight(0.6));
    EXPECT_DOUBLE_EQ(data_d0.correlated_probabilities_sum, 0.6);
    ASSERT_EQ(data_d0.implied_weights_for_other_edges.size(), 1);
    pm::ImpliedWeightUnconverted expected_iwu_d0_to_d1 = {d1.d1, d1.d2, pm::convert_probability_to_weight(0.9)};
    EXPECT_EQ(data_d0.implied_weights_for_other_edges[0], expected_iwu_d0_to_d1);

    // Check D1
    ASSERT_TRUE(graph.edge_map.count(d1));
    const auto& data_d1 = graph.edge_map.at(d1);
    EXPECT_EQ(data_d1.detectors, d1);
    EXPECT_EQ(data_d1.observables, std::vector<size_t>{0});
    EXPECT_DOUBLE_EQ(data_d1.error_probability, 0.6);
    EXPECT_DOUBLE_EQ(data_d1.weight, convert_p_to_double_weight(0.6));
    ASSERT_EQ(data_d1.implied_weights_for_other_edges.size(), 1);
    pm::ImpliedWeightUnconverted expected_iwu_d1_to_d0 = {d0.d1, d0.d2, pm::convert_probability_to_weight(0.9)};
    EXPECT_EQ(data_d1.implied_weights_for_other_edges[0], expected_iwu_d1_to_d0);
}

TEST(DetectorErrorModelToUserGraph, ErrorWithNoDetectorsInComponent) {
    stim::DetectorErrorModel dem("error(0.1) L0 ^ D1");
    pm::UserGraph graph = pm::detector_error_model_to_user_graph(dem);

    pm::DetectorEdgeId d1(1);
    ASSERT_EQ(graph.edge_map.size(), 1);  // Only D1 edge

    // Check D1
    ASSERT_TRUE(graph.edge_map.count(d1));
    const auto& data_d1 = graph.edge_map.at(d1);
    EXPECT_EQ(data_d1.detectors, d1);
    EXPECT_TRUE(data_d1.observables.empty());
    EXPECT_DOUBLE_EQ(data_d1.error_probability, 0.1);
    EXPECT_DOUBLE_EQ(data_d1.weight, convert_p_to_double_weight(0.1));
    EXPECT_DOUBLE_EQ(data_d1.correlated_probabilities_sum, 0.1);
    EXPECT_TRUE(data_d1.implied_weights_for_other_edges.empty());
}

TEST(DetectorErrorModelToUserGraph, ComplexInteractingCorrelatedErrors) {
    const char* dem_text = R"DEM(
        error(0.1) L0 D0 D1 ^ D2        #I1
        error(0.2) D0 D1 ^ D2 ^ D3   #I2
        error(0.15) D3 ^ D0 D1       #I3
    )DEM";
    stim::DetectorErrorModel dem(dem_text);
    pm::UserGraph graph = pm::detector_error_model_to_user_graph(dem);

    pm::DetectorEdgeId d01(0, 1);
    pm::DetectorEdgeId d2(2);
    pm::DetectorEdgeId d3(3);

    ASSERT_EQ(graph.edge_map.size(), 3);  // D(0,1), D(2), D(3)

    // --- Expected error_probability and weight (from handle_dem_instruction calls) ---
    // D(0,1): Calls with p=0.1 (I1), p=0.2 (I2), p=0.15 (I3)
    // ep_d01_i1 = 0.1
    // ep_d01_i2 = bernoulli_xor(0.1, 0.2) = 0.26
    // ep_d01_i3 = bernoulli_xor(0.26, 0.15) = 0.26*0.85 + 0.15*0.74 = 0.221 + 0.111 = 0.332
    double ep_d01 = 0.332;

    // D(2): Calls with p=0.1 (I1), p=0.2 (I2)
    // ep_d2_i1 = 0.1
    // ep_d2_i2 = bernoulli_xor(0.1, 0.2) = 0.26
    double ep_d2 = 0.26;

    // D(3): Calls with p=0.2 (I2), p=0.15 (I3)
    // ep_d3_i1 = 0.2 // (from I2)
    // ep_d3_i2 = bernoulli_xor(0.2, 0.15) = 0.2*0.85 + 0.15*0.8 = 0.17 + 0.12 = 0.29
    double ep_d3 = 0.29;

    // --- Expected conditional_groups ---
    // cg = { d01:{d2:0.26, d3:0.29}, d2:{d01:0.26, d3:0.2}, d3:{d01:0.29, d2:0.2} }

    // --- Check D(0,1) ---
    ASSERT_TRUE(graph.edge_map.count(d01));
    const auto& data_d01 = graph.edge_map.at(d01);
    EXPECT_EQ(data_d01.detectors, d01);
    EXPECT_EQ(data_d01.observables.size(), 1);
    EXPECT_EQ(data_d01.observables[0], 0);
    EXPECT_DOUBLE_EQ(data_d01.error_probability, ep_d01);
    EXPECT_DOUBLE_EQ(data_d01.weight, convert_p_to_double_weight(ep_d01));
    // Causal d01: sum_prob = cg[d01][d2] + cg[d01][d3] = 0.26 + 0.29 = 0.55
    EXPECT_DOUBLE_EQ(data_d01.correlated_probabilities_sum, 0.55);
    ASSERT_EQ(data_d01.implied_weights_for_other_edges.size(), 2);
    // d01->d2: impl_p=min(0.9, 0.26/0.55) approx 0.472727
    // d01->d3: impl_p=min(0.9, 0.29/0.55) approx 0.527272
    pm::ImpliedWeightUnconverted expected_d01_to_d2 = {d2.d1, d2.d2, pm::convert_probability_to_weight(0.26 / 0.55)};
    pm::ImpliedWeightUnconverted expected_d01_to_d3 = {d3.d1, d3.d2, pm::convert_probability_to_weight(0.29 / 0.55)};
    EXPECT_EQ(data_d01.implied_weights_for_other_edges[0], expected_d01_to_d2);  // d2 < d3
    EXPECT_EQ(data_d01.implied_weights_for_other_edges[1], expected_d01_to_d3);

    // --- Check D(2) ---
    ASSERT_TRUE(graph.edge_map.count(d2));
    const auto& data_d2 = graph.edge_map.at(d2);
    EXPECT_EQ(data_d2.detectors, d2);
    EXPECT_TRUE(data_d2.observables.empty());
    EXPECT_DOUBLE_EQ(data_d2.error_probability, ep_d2);
    EXPECT_DOUBLE_EQ(data_d2.weight, convert_p_to_double_weight(ep_d2));
    // Causal d2: sum_prob = cg[d2][d01] + cg[d2][d3] = 0.26 + 0.2 = 0.46
    EXPECT_DOUBLE_EQ(data_d2.correlated_probabilities_sum, 0.46);
    ASSERT_EQ(data_d2.implied_weights_for_other_edges.size(), 2);
    // d2->d01: impl_p=min(0.9, 0.26/0.46) approx 0.565217
    // d2->d3:  impl_p=min(0.9, 0.2/0.46) approx 0.434782
    pm::ImpliedWeightUnconverted expected_d2_to_d01 = {d01.d1, d01.d2, pm::convert_probability_to_weight(0.26 / 0.46)};
    pm::ImpliedWeightUnconverted expected_d2_to_d3 = {d3.d1, d3.d2, pm::convert_probability_to_weight(0.2 / 0.46)};
    EXPECT_EQ(data_d2.implied_weights_for_other_edges[0], expected_d2_to_d01);  // d01 < d3
    EXPECT_EQ(data_d2.implied_weights_for_other_edges[1], expected_d2_to_d3);

    // --- Check D(3) ---
    ASSERT_TRUE(graph.edge_map.count(d3));
    const auto& data_d3 = graph.edge_map.at(d3);
    EXPECT_EQ(data_d3.detectors, d3);
    EXPECT_TRUE(data_d3.observables.empty());
    EXPECT_DOUBLE_EQ(data_d3.error_probability, ep_d3);
    EXPECT_DOUBLE_EQ(data_d3.weight, convert_p_to_double_weight(ep_d3));
    // Causal d3: sum_prob = cg[d3][d01] + cg[d3][d2] = 0.29 + 0.2 = 0.49
    EXPECT_DOUBLE_EQ(data_d3.correlated_probabilities_sum, 0.49);
    ASSERT_EQ(data_d3.implied_weights_for_other_edges.size(), 2);
    // d3->d01: impl_p=min(0.9, 0.29/0.49) approx 0.591836
    // d3->d2:  impl_p=min(0.9, 0.2/0.49) approx 0.408163
    pm::ImpliedWeightUnconverted expected_d3_to_d01 = {d01.d1, d01.d2, pm::convert_probability_to_weight(0.29 / 0.49)};
    pm::ImpliedWeightUnconverted expected_d3_to_d2 = {d2.d1, d2.d2, pm::convert_probability_to_weight(0.2 / 0.49)};
    EXPECT_EQ(data_d3.implied_weights_for_other_edges[0], expected_d3_to_d01);  // d01 < d2
    EXPECT_EQ(data_d3.implied_weights_for_other_edges[1], expected_d3_to_d2);
}
