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

#include "pymatching/sparse_blossom/driver/io.h"

#include "gtest/gtest.h"

using pm::DetectorEdgeId;
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

/// Tests for data structures used within the UserGraph.
TEST(DetectorEdgeIdOperatorsTest, Construction) {
    DetectorEdgeId dempty;
    EXPECT_EQ(dempty.d1, SIZE_MAX);
    EXPECT_EQ(dempty.d2, SIZE_MAX);

    DetectorEdgeId d0(0);
    EXPECT_EQ(d0.d1, 0);
    EXPECT_EQ(d0.d2, SIZE_MAX);

    DetectorEdgeId d12(1, 2);
    EXPECT_EQ(d12.d1, 1);
    EXPECT_EQ(d12.d2, 2);
}

TEST(DetectorEdgeIdOperatorsTest, Equality) {
    DetectorEdgeId d12(1, 2);
    DetectorEdgeId d21(1, 2);
    DetectorEdgeId d22(2, 2);
    EXPECT_EQ(d12, d21);
    EXPECT_NE(d22, d21);
    EXPECT_NE(d22, d12);
}

TEST(DetectorEdgeIdOperatorsTest, IsValidEdge) {
    DetectorEdgeId dinvalid;
    DetectorEdgeId d0(0);
    DetectorEdgeId d12(1, 2);
    EXPECT_FALSE(dinvalid.is_valid_edge());
    EXPECT_TRUE(d0.is_valid_edge());
    EXPECT_TRUE(d12.is_valid_edge());
}

TEST(DetectorEdgeIdOperatorsTest, IsBoundary) {
    DetectorEdgeId dinvalid;
    DetectorEdgeId d0(0);
    DetectorEdgeId d12(1, 2);
    EXPECT_FALSE(dinvalid.is_boundary());
    EXPECT_TRUE(d0.is_boundary());
    EXPECT_FALSE(d12.is_boundary());
}

TEST(DetectorEdgeDataOperatorsTest, EqualityAndInequality) {
    pm::DetectorEdgeData d1{1, {10, 20}, 0.5, 0.01};
    pm::DetectorEdgeData d2{1, {10, 20}, 0.5, 0.01};  // Equal to d1
    pm::DetectorEdgeData d3{2, {10, 20}, 0.5, 0.01};  // Different detectors
    pm::DetectorEdgeData d4{1, {10, 21}, 0.5, 0.01};  // Different observables
    pm::DetectorEdgeData d5{1, {10, 20}, 0.6, 0.01};  // Different weight
    pm::DetectorEdgeData d6{1, {10, 20}, 0.5, 0.02};  // Different error_probability

    EXPECT_TRUE(d1 == d2);
    EXPECT_FALSE(d1 != d2);

    EXPECT_FALSE(d1 == d3);
    EXPECT_TRUE(d1 != d3);

    EXPECT_FALSE(d1 == d4);
    EXPECT_TRUE(d1 != d4);

    EXPECT_FALSE(d1 == d5);
    EXPECT_TRUE(d1 != d5);

    EXPECT_FALSE(d1 == d6);
    EXPECT_TRUE(d1 != d6);

    // Test when a vector member is empty vs non-empty but otherwise identical fields
    pm::DetectorEdgeData d7{1, {}, 0.5, 0.01};
    pm::DetectorEdgeData d8{1, {10}, 0.5, 0.01};
    EXPECT_FALSE(d7 == d1);  // d1 has observables {10,20}
    EXPECT_TRUE(d7 != d1);
    EXPECT_FALSE(d7 == d8);
    EXPECT_TRUE(d7 != d8);
}

/// Tests for DEM iteration.

struct HandlerCall {
    double p;
    DetectorEdgeId dets;
    std::vector<size_t> observables;

    HandlerCall(double probability, DetectorEdgeId detectors, std::vector<size_t> obs)
        : p(probability), dets(detectors), observables(std::move(obs)) {
    }

    bool operator==(const HandlerCall& other) const {
        auto sorted_obs_this = this->observables;
        auto sorted_obs_other = other.observables;
        std::sort(sorted_obs_this.begin(), sorted_obs_this.end());
        std::sort(sorted_obs_other.begin(), sorted_obs_other.end());

        return p == other.p && dets == other.dets && sorted_obs_this == sorted_obs_other;
    }

    friend std::ostream& operator<<(std::ostream& os, const HandlerCall& call) {
        os << "HandlerCall{p=" << call.p << ", dets=" << call.dets.d1 << ", " << call.dets.d2 << ", obs=[";
        for (size_t i = 0; i < call.observables.size(); ++i) {
            os << call.observables[i] << (i == call.observables.size() - 1 ? "" : ", ");
        }
        os << "]}";
        return os;
    }
};

TEST(IterDemInstructionsIncludeCorrelations, EmptyDem) {
    stim::DetectorErrorModel dem;  // Empty model
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    EXPECT_TRUE(handler_calls.empty());
    EXPECT_TRUE(conditional_groups.empty());
}

TEST(IterDemInstructionsIncludeCorrelations, SingleDetectorError) {
    const char* single_det_dem = "error(0.1) D0";
    stim::DetectorErrorModel dem(single_det_dem);
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> expected_calls = {{0.1, DetectorEdgeId(0), {}}};
    std::vector<HandlerCall> handler_calls;

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), 1);
    ASSERT_EQ(handler_calls, expected_calls);

    // Validate conditional_groups
    // Decomposed error: prob=0.1, components=[{D(0)}]
    // Result: cg[D(0)][D(0)] = 0.1
    DetectorEdgeId d0(0);
    EXPECT_EQ(conditional_groups.size(), 1);
    EXPECT_TRUE(conditional_groups.count(d0));
    EXPECT_EQ(conditional_groups.at(d0).size(), 1);
    EXPECT_EQ(conditional_groups.at(d0).count(d0), 1);
    EXPECT_DOUBLE_EQ(conditional_groups.at(d0).at(d0), 0.1);
}

TEST(IterDemInstructionsIncludeCorrelations, TwoDetectorError) {
    stim::DetectorErrorModel dem("error(0.2) D1 D3");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> expected_calls = {{0.2, DetectorEdgeId(1, 3), {}}};
    std::vector<HandlerCall> handler_calls;

    // Handler signature uses std::vector<size_t>
    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), 1);
    ASSERT_EQ(handler_calls, expected_calls);

    // Validate conditional_groups
    // Decomposed error: prob=0.2, components=[{dets=D(1,3)}]
    // Result: cg[D(1,3)][D(1,3)] = 0.2
    DetectorEdgeId d13(1, 3);
    EXPECT_EQ(conditional_groups.size(), 1);
    ASSERT_TRUE(conditional_groups.count(d13));
    EXPECT_EQ(conditional_groups.at(d13).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d13).count(d13));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d13).at(d13), 0.2);
}

TEST(IterDemInstructionsIncludeCorrelations, ThreeDetectorErrorIsIgnoredByHandler) {
    stim::DetectorErrorModel dem("error(0.3) D0 D1 D2");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    EXPECT_TRUE(handler_calls.empty());
    EXPECT_TRUE(conditional_groups.empty());
}

TEST(IterDemInstructionsIncludeCorrelations, ThreeDetectorErrorWithObservableIsIgnoredByHandler) {
    stim::DetectorErrorModel dem("error(0.3) L0 D0 D1 D2");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    EXPECT_TRUE(handler_calls.empty());
    ASSERT_TRUE(conditional_groups.empty());
}

TEST(IterDemInstructionsIncludeCorrelations, SingleDetectorWithObservable) {
    stim::DetectorErrorModel dem("error(0.4) D5 L2");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<HandlerCall> expected_calls = {{0.4, DetectorEdgeId(5), {2}}};

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), 1);
    ASSERT_EQ(handler_calls, expected_calls);

    // Validate conditional_groups
    // Decomposed error: prob=0.4, components=[{dets=D(5), obs={2}}]
    // Result: cg[D(5)][D(5)] = 0.4
    DetectorEdgeId d5(5);
    EXPECT_EQ(conditional_groups.size(), 1);
    ASSERT_TRUE(conditional_groups.count(d5));
    EXPECT_EQ(conditional_groups.at(d5).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d5).count(d5));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d5).at(d5), 0.4);
}

TEST(IterDemInstructionsIncludeCorrelations, TwoDetectorsWithObservables) {
    stim::DetectorErrorModel dem("error(0.5) D1 D3 L0 L4");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<HandlerCall> expected_calls = {{0.5, DetectorEdgeId(1, 3), {0, 4}}};

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), 1);
    ASSERT_EQ(handler_calls, expected_calls);
    // Validate conditional_groups
    // Decomposed error: prob=0.5, components=[{dets=D(1,3), obs={0,4}}]
    // Result: cg[D(1,3)][D(1,3)] = 0.5
    DetectorEdgeId d13(1, 3);
    EXPECT_EQ(conditional_groups.size(), 1);
    ASSERT_TRUE(conditional_groups.count(d13));
    EXPECT_EQ(conditional_groups.at(d13).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d13).count(d13));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d13).at(d13), 0.5);
}

TEST(IterDemInstructionsIncludeCorrelations, ZeroProbabilityError) {
    stim::DetectorErrorModel dem("error(0.0) D0");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    EXPECT_TRUE(handler_calls.empty());
}

TEST(IterDemInstructionsIncludeCorrelations, CorrelatedErrorSimple) {
    stim::DetectorErrorModel dem("error(0.6) D0 ^ D1 L0");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<HandlerCall> expected_calls = {{0.6, DetectorEdgeId(0), {}}, {0.6, DetectorEdgeId(1), {0}}};

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), expected_calls.size());
    EXPECT_EQ(handler_calls, expected_calls);
    // Validate conditional_groups
    // Decomposed error: prob=0.6, components=[{dets=D(0)}, {dets=D(1), obs={0}}]
    // Result: cg[D(0)][D(1)] = 0.6, cg[D(1)][D(0)] = 0.6
    DetectorEdgeId d0(0), d1(1);
    EXPECT_EQ(conditional_groups.size(), 2);
    ASSERT_TRUE(conditional_groups.count(d0));
    ASSERT_TRUE(conditional_groups.count(d1));

    EXPECT_EQ(conditional_groups.at(d0).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d0).count(d1));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d0).at(d1), 0.6);

    EXPECT_EQ(conditional_groups.at(d1).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d1).count(d0));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d1).at(d0), 0.6);
}

TEST(IterDemInstructionsIncludeCorrelations, CorrelatedErrorThreeComponents) {
    stim::DetectorErrorModel dem("error(0.7) D0 L1 ^ D2 D3 L0 ^ D4");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<HandlerCall> expected_calls = {
        {0.7, DetectorEdgeId(0), {1}}, {0.7, DetectorEdgeId(2, 3), {0}}, {0.7, DetectorEdgeId(4), {}}};

    // Handler signature uses std::vector<size_t>
    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), expected_calls.size());
    EXPECT_EQ(handler_calls, expected_calls);

    // Validate conditional_groups
    // Decomposed error: prob=0.7, components=[{dets=D(0),obs={1}}, {dets=D(2,3),obs={0}}, {dets=D(4)}]
    // Pairs: (D0,D23), (D0,D4), (D23,D4), all get probability 0.7
    DetectorEdgeId d0(0), d23(2, 3), d4(4);
    EXPECT_EQ(conditional_groups.size(), 3);
    ASSERT_TRUE(conditional_groups.count(d0));
    ASSERT_TRUE(conditional_groups.count(d23));
    ASSERT_TRUE(conditional_groups.count(d4));

    // D0 pairs
    EXPECT_EQ(conditional_groups.at(d0).size(), 2);
    ASSERT_TRUE(conditional_groups.at(d0).count(d23));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d0).at(d23), 0.7);
    ASSERT_TRUE(conditional_groups.at(d0).count(d4));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d0).at(d4), 0.7);

    // D23 pairs
    EXPECT_EQ(conditional_groups.at(d23).size(), 2);
    ASSERT_TRUE(conditional_groups.at(d23).count(d0));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d23).at(d0), 0.7);
    ASSERT_TRUE(conditional_groups.at(d23).count(d4));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d23).at(d4), 0.7);

    // D4 pairs
    EXPECT_EQ(conditional_groups.at(d4).size(), 2);
    ASSERT_TRUE(conditional_groups.at(d4).count(d0));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d4).at(d0), 0.7);
    ASSERT_TRUE(conditional_groups.at(d4).count(d23));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d4).at(d23), 0.7);
}

TEST(IterDemInstructionsIncludeCorrelations, CorrelatedErrorWithBoundaryComponentLast) {
    stim::DetectorErrorModel dem("error(0.8) D0 D1 D2 ^ D3 L1");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<HandlerCall> expected_calls = {{0.8, DetectorEdgeId(3), {1}}};

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), expected_calls.size());
    EXPECT_EQ(handler_calls, expected_calls);
    // Validate conditional_groups
    // Boundary D0D1D2 removed. Decomposed: prob=0.8, components=[{dets=D(3), obs={1}}]
    // Result: cg[D(3)][D(3)] = 0.8
    DetectorEdgeId d3(3);
    EXPECT_EQ(conditional_groups.size(), 1);
    ASSERT_TRUE(conditional_groups.count(d3));
    EXPECT_EQ(conditional_groups.at(d3).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d3).count(d3));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d3).at(d3), 0.8);
}

TEST(IterDemInstructionsIncludeCorrelations, CorrelatedErrorWithBoundaryComponentFirst) {
    stim::DetectorErrorModel dem("error(0.9) D0 L0 ^ D1 D2 D3");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<HandlerCall> expected_calls = {{0.9, DetectorEdgeId(0), {0}}};

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), expected_calls.size());
    EXPECT_EQ(handler_calls, expected_calls);

    // Validate conditional_groups
    // Boundary D1D2D3 removed. Decomposed: prob=0.9, components=[{dets=D(0), obs={0}}]
    // Result: cg[D(0)][D(0)] = 0.9
    DetectorEdgeId d0(0);
    EXPECT_EQ(conditional_groups.size(), 1);
    ASSERT_TRUE(conditional_groups.count(d0));
    EXPECT_EQ(conditional_groups.at(d0).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d0).count(d0));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d0).at(d0), 0.9);
}

TEST(IterDemInstructionsIncludeCorrelations, HandlesDecompositionElemWithNoDetectors) {
    stim::DetectorErrorModel dem("error(0.1) L0 ^ D1");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<HandlerCall> expected_calls = {{0.1, DetectorEdgeId(1), {}}};

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);
    ASSERT_EQ(handler_calls.size(), expected_calls.size());
    EXPECT_EQ(handler_calls, expected_calls);
    // Validate conditional_groups
    // First component (L0 with no detectors) is invalid and popped.
    // Decomposed error: prob=0.1, components=[{dets=D(1)}]
    // Result: cg[D(1)][D(1)] = 0.1
    DetectorEdgeId d1(1);
    EXPECT_EQ(conditional_groups.size(), 1);
    ASSERT_TRUE(conditional_groups.count(d1));
    EXPECT_EQ(conditional_groups.at(d1).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d1).count(d1));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d1).at(d1), 0.1);
}

TEST(IterDemInstructionsIncludeCorrelations, HandlesDecompositionElemWithNoDetectorsAfterSeparator) {
    stim::DetectorErrorModel dem("error(0.1) D0 ^ L1");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<HandlerCall> expected_calls = {{0.1, DetectorEdgeId(0), {}}};

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);
    ASSERT_EQ(handler_calls.size(), expected_calls.size());
    EXPECT_EQ(handler_calls, expected_calls);
    EXPECT_EQ(conditional_groups.size(), 1);

    // Validate conditional_groups
    // Second component (L1 with no detectors) is invalid and popped.
    // Decomposed error: prob=0.1, components=[{dets=D(0)}]
    // Result: cg[D(0)][D(0)] = 0.1
    DetectorEdgeId d0(0);
    EXPECT_EQ(conditional_groups.size(), 1);
    ASSERT_TRUE(conditional_groups.count(d0));
    EXPECT_EQ(conditional_groups.at(d0).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d0).count(d0));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d0).at(d0), 0.1);
}

TEST(IterDemInstructionsIncludeCorrelations, MultipleInstructions) {
    const char* dem_text = R"DEM(
        error(0.1) D0           # Instruction 1: Simple
        error(0.2) D1 D2 L0    # Instruction 2: Two detectors, one observable
        error(0.3) D3 D4 D5 ^ D6 # Instruction 3: Boundary ignored, second component handled
        error(0.0) D7           # Instruction 4: Zero probability, ignored
        error(0.4) D8 ^ D9 L1    # Instruction 5: Correlated
    )DEM";
    stim::DetectorErrorModel dem(dem_text);
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<HandlerCall> expected_calls = {
        // From instruction 1
        {0.1, DetectorEdgeId(0), {}},  // obs {}
        // From instruction 2
        {0.2, DetectorEdgeId(1, 2), {0}},  // obs {0}
        // From instruction 3 (only D6 part)
        {0.3, DetectorEdgeId(6), {}},  // obs {}
        // From instruction 4 (none, p=0)
        // From instruction 5
        {0.4, DetectorEdgeId(8), {}},   // obs {}
        {0.4, DetectorEdgeId(9), {1}},  // obs {1}
    };

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), expected_calls.size());
    EXPECT_EQ(handler_calls, expected_calls);

    // Validate conditional_groups
    DetectorEdgeId d0(0), d12(1, 2), d6(6), d7(7), d8(8), d9(9);

    // Expected states after each instruction processing by add_decomposed_error_to_conditional_groups:
    // I1 (error(0.1) D0): dec_err={prob=0.1, comp=[{D0}]}. cg[D0][D0] = 0.1
    // I2 (error(0.2) D1 D2 L0): dec_err={prob=0.2, comp=[{D12}]}. cg[D12][D12] = 0.2
    // I3 (error(0.3) D3 D4 D5 ^ D6): dec_err={prob=0.3, comp=[{D6}]}. cg[D6][D6] = 0.3
    // I4 (error(0.0) D7): dec_err={prob=0.0, comp=[{D7}]}. cg[D7][D7] = 0.0
    // I5 (error(0.4) D8 ^ D9 L1): dec_err={prob=0.4, comp=[{D8},{D9,L1}]}. cg[D8][D9]=0.4, cg[D9][D8]=0.4

    EXPECT_EQ(conditional_groups.size(), 6);  // Keys: D0, D(1,2), D6, D7, D8, D9

    // Check D0
    ASSERT_TRUE(conditional_groups.count(d0));
    EXPECT_EQ(conditional_groups.at(d0).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d0).count(d0));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d0).at(d0), 0.1);

    // Check D(1,2)
    ASSERT_TRUE(conditional_groups.count(d12));
    EXPECT_EQ(conditional_groups.at(d12).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d12).count(d12));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d12).at(d12), 0.2);

    // Check D6
    ASSERT_TRUE(conditional_groups.count(d6));
    EXPECT_EQ(conditional_groups.at(d6).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d6).count(d6));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d6).at(d6), 0.3);

    // Check D7
    ASSERT_TRUE(conditional_groups.count(d7));
    EXPECT_EQ(conditional_groups.at(d7).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d7).count(d7));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d7).at(d7), 0.0);

    // Check D8 and D9 (coupled from I5)
    ASSERT_TRUE(conditional_groups.count(d8));
    EXPECT_EQ(conditional_groups.at(d8).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d8).count(d9));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d8).at(d9), 0.4);

    ASSERT_TRUE(conditional_groups.count(d9));
    EXPECT_EQ(conditional_groups.at(d9).size(), 1);
    ASSERT_TRUE(conditional_groups.at(d9).count(d8));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d9).at(d8), 0.4);
}

TEST(IterDemInstructionsIncludeCorrelations, ComplexInteractingCorrelatedErrors) {
    const char* dem_text = R"DEM(
        error(0.1) D0 D1 ^ D2
        error(0.2) D0 D1 ^ D2 ^ D3
        error(0.15) D3 ^ D0 D1
    )DEM";
    stim::DetectorErrorModel dem(dem_text);
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;

    DetectorEdgeId d01(0, 1);
    DetectorEdgeId d2(2);
    DetectorEdgeId d3(3);

    std::vector<HandlerCall> expected_calls = {
        // From error(0.1) D0 D1 ^ D2
        {0.1, d01, {}},
        {0.1, d2, {}},
        // From error(0.2) D0 D1 ^ D2 ^ D3
        {0.2, d01, {}},
        {0.2, d2, {}},
        {0.2, d3, {}},
        // From error(0.15) D3 ^ D0 D1
        {0.15, d3, {}},
        {0.15, d01, {}}};

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), expected_calls.size());
    EXPECT_EQ(handler_calls, expected_calls);

    // Validate conditional_groups step-by-step calculation:
    // Initial: cg = {}
    //
    // After I1 (error(0.1) D0 D1 ^ D2):
    //   decomposed_err = {prob=0.1, comp=[{d01}, {d2}]}
    //   cg[d01][d2] = bernoulli_xor(0.0, 0.1) = 0.1
    //   cg[d2][d01] = bernoulli_xor(0.0, 0.1) = 0.1
    //   Current cg: { d01:{d2:0.1}, d2:{d01:0.1} }
    //
    // After I2 (error(0.2) D0 D1 ^ D2 ^ D3):
    //   decomposed_err = {prob=0.2, comp=[{d01}, {d2}, {d3}]}
    //   Pairs: (d01,d2), (d01,d3), (d2,d3) with p=0.2
    //   cg[d01][d2] = bernoulli_xor(0.1, 0.2) = 0.1*0.8 + 0.2*0.9 = 0.08 + 0.18 = 0.26
    //   cg[d2][d01] = bernoulli_xor(0.1, 0.2) = 0.26
    //   cg[d01][d3] = bernoulli_xor(0.0, 0.2) = 0.2
    //   cg[d3][d01] = bernoulli_xor(0.0, 0.2) = 0.2
    //   cg[d2][d3]  = bernoulli_xor(0.0, 0.2) = 0.2
    //   cg[d3][d2]  = bernoulli_xor(0.0, 0.2) = 0.2
    //   Current cg: { d01:{d2:0.26, d3:0.2}, d2:{d01:0.26, d3:0.2}, d3:{d01:0.2, d2:0.2} }
    //
    // After I3 (error(0.15) D3 ^ D0 D1):
    //   decomposed_err = {prob=0.15, comp=[{d3}, {d01}]}
    //   Pair: (d3,d01) with p=0.15
    //   cg[d3][d01]  = bernoulli_xor(0.2, 0.15) = 0.2*0.85 + 0.15*0.8 = 0.17 + 0.12 = 0.29
    //   cg[d01][d3]  = bernoulli_xor(0.2, 0.15) = 0.29
    //
    // Final expected state for conditional_groups:
    // {
    //   d01: { d2: 0.26, d3: 0.29 },
    //   d2:  { d01: 0.26, d3: 0.2 },
    //   d3:  { d01: 0.29, d2: 0.2 }
    // }

    EXPECT_EQ(conditional_groups.size(), 3);  // Keys: d01, d2, d3

    // Check for d01
    ASSERT_TRUE(conditional_groups.count(d01));
    EXPECT_EQ(conditional_groups.at(d01).size(), 2);
    ASSERT_TRUE(conditional_groups.at(d01).count(d2));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d01).at(d2), 0.26);
    ASSERT_TRUE(conditional_groups.at(d01).count(d3));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d01).at(d3), 0.29);

    // Check for d2
    ASSERT_TRUE(conditional_groups.count(d2));
    EXPECT_EQ(conditional_groups.at(d2).size(), 2);
    ASSERT_TRUE(conditional_groups.at(d2).count(d01));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d2).at(d01), 0.26);
    ASSERT_TRUE(conditional_groups.at(d2).count(d3));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d2).at(d3), 0.2);

    // Check for d3
    ASSERT_TRUE(conditional_groups.count(d3));
    EXPECT_EQ(conditional_groups.at(d3).size(), 2);
    ASSERT_TRUE(conditional_groups.at(d3).count(d01));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d3).at(d01), 0.29);
    ASSERT_TRUE(conditional_groups.at(d3).count(d2));
    EXPECT_DOUBLE_EQ(conditional_groups.at(d3).at(d2), 0.2);
}
