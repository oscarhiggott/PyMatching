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
    std::vector<HandlerCall> handler_calls;

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), 1);

    EXPECT_DOUBLE_EQ(handler_calls[0].p, 0.1);
    EXPECT_EQ(handler_calls[0].dets, DetectorEdgeId(0));
    EXPECT_TRUE(handler_calls[0].observables.empty());
}

TEST(IterDemInstructionsIncludeCorrelations, TwoDetectorError) {
    stim::DetectorErrorModel dem("error(0.2) D1 D3");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;

    // Handler signature uses std::vector<size_t>
    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), 1);
    EXPECT_DOUBLE_EQ(handler_calls[0].p, 0.2);
    EXPECT_EQ(handler_calls[0].dets, DetectorEdgeId(1, 3));
    EXPECT_TRUE(handler_calls[0].observables.empty());
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
}

TEST(IterDemInstructionsIncludeCorrelations, ThreeDetectorErrorWithObservableIsIgnoredByHandler) {
    stim::DetectorErrorModel dem("error(0.3) D0 D1 D2");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    EXPECT_TRUE(handler_calls.empty());
}

TEST(IterDemInstructionsIncludeCorrelations, SingleDetectorWithObservable) {
    stim::DetectorErrorModel dem("error(0.4) D5 L2");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<size_t> expected_obs = {2};

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), 1);
    EXPECT_DOUBLE_EQ(handler_calls[0].p, 0.4);
    EXPECT_EQ(handler_calls[0].dets, DetectorEdgeId(5));
    EXPECT_EQ(handler_calls[0].observables, expected_obs);
}

TEST(IterDemInstructionsIncludeCorrelations, TwoDetectorsWithObservables) {
    stim::DetectorErrorModel dem("error(0.5) D1 D3 L0 L4");
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> conditional_groups;
    std::vector<HandlerCall> handler_calls;
    std::vector<size_t> expected_obs = {0, 4};

    auto handler = [&](double p, const DetectorEdgeId& dets, const std::vector<size_t>& obs) {
        handler_calls.emplace_back(p, dets, obs);
    };

    iter_dem_instructions_include_correlations(dem, handler, conditional_groups);

    ASSERT_EQ(handler_calls.size(), 1);
    EXPECT_DOUBLE_EQ(handler_calls[0].p, 0.5);
    EXPECT_EQ(handler_calls[0].dets, DetectorEdgeId(1, 3));
    EXPECT_EQ(handler_calls[0].observables, expected_obs);
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
}
