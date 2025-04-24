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

#ifndef PYMATCHING2_IO_H
#define PYMATCHING2_IO_H

#include "pymatching/sparse_blossom/flooder/graph.h"
#include "stim.h"

namespace pm {

const pm::weight_int NUM_DISTINCT_WEIGHTS = 1 << (sizeof(pm::weight_int) * 8 - 8);

/// Computes the weight of an edge resulting from merging edges with weight `a' and weight `b', assuming each edge
/// weight is a log-likelihood ratio log((1-p)/p) associated with the probability p of an error occurring on the
/// edge, and that the error mechanisms associated with the two edges being merged are independent.
///
/// A mathematically equivalent implementation of this method would be:
///
///  double merge_weights(double a, double b){
///     double p_a = 1/(1 + std::exp(a));
///     double p_b = 1/(1 + std::exp(b));
///     double p_both = p_a * (1 - p_b) + p_b * (1 - p_a);
///     return std::log((1-p_both)/p_both);
///  }
///
/// however this would suffer from numerical overflow issues for abs(a) >> 30 or abs(b) >> 30 which is avoided by the
/// the implementation used here instead. See Equation (6) of https://ieeexplore.ieee.org/document/1495850 for more
/// details.
double merge_weights(double a, double b);

// Identifies an edge defined by two detectors (d1, d2) in a detector graph.
// d1 <= d2 for convenience is guaranteed for convenience.
//
// If d2 is `SIZE_MAX` and d1 < SIZE_MAX, then the edge ID is considered to be that of a boundary edge.
//
// An edge is considered a valid edge if it has at least one detector assigned to it.
struct DetectorEdgeId {
    size_t d1;
    size_t d2;

    DetectorEdgeId();
    DetectorEdgeId(const size_t d);
    DetectorEdgeId(const size_t d1, const size_t d2);

    bool operator==(const DetectorEdgeId &other) const;
    bool operator!=(const DetectorEdgeId &other) const;
    bool operator<(const DetectorEdgeId &other) const;

    bool is_valid_edge() const;
    bool is_boundary() const;
};

struct DetectorEdgeData {
    DetectorEdgeId detectors;
    std::vector<size_t> observables;      /// The observables crossed along this edge
    double weight;                        /// The weight of the edge to this neighboring node
    double error_probability;             /// The error probability associated with this node
    double correlated_probabilities_sum;  /// The error probability associated with this node
    std::vector<ImpliedWeightUnconverted> implied_weights_for_other_edges{};
    bool operator==(const DetectorEdgeData &other) const;
    bool operator!=(const DetectorEdgeData &other) const;
};

struct DecomposedDemError {
    /// The probability of this error occurring.
    double probability;
    /// Effects of the error.
    stim::FixedCapVector<DetectorEdgeData, 8> components;

    bool operator==(const DecomposedDemError &other) const;
    bool operator!=(const DecomposedDemError &other) const;
};

void add_decomposed_error_to_conditional_groups(
    pm::DecomposedDemError &error, std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> &conditional_groups);

template <typename Handler>
void iter_detector_error_model_edges(
    const stim::DetectorErrorModel &detector_error_model, const Handler &handle_dem_error) {
    detector_error_model.iter_flatten_error_instructions([&](const stim::DemInstruction &instruction) {
        std::vector<size_t> dets;
        std::vector<size_t> observables;
        double p = instruction.arg_data[0];
        for (auto &target : instruction.target_data) {
            if (target.is_relative_detector_id()) {
                dets.push_back(target.val());
            } else if (target.is_observable_id()) {
                observables.push_back(target.val());
            } else if (target.is_separator()) {
                if (p > 0) {
                    handle_dem_error(p, dets, observables);
                    observables.clear();
                    dets.clear();
                }
            }
        }
        if (p > 0) {
            handle_dem_error(p, dets, observables);
        }
    });
}

template <typename Handler>
void iter_dem_instructions_include_correlations(
    const stim::DetectorErrorModel &detector_error_model,
    const Handler &handle_dem_error,
    std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> &conditional_groups) {
    detector_error_model.iter_flatten_error_instructions([&](const stim::DemInstruction &instruction) {
        double p = instruction.arg_data[0];
        pm::DecomposedDemError decomposed_err;
        decomposed_err.probability = p;
        decomposed_err.components = {};
        decomposed_err.components.push_back({});
        DetectorEdgeData *component = &decomposed_err.components.back();
        size_t num_translated_relative_detectors = 0;
        for (auto &target : instruction.target_data) {
            // Decompose error
            if (target.is_relative_detector_id()) {
                num_translated_relative_detectors++;
                if (num_translated_relative_detectors == 1) {
                    const size_t &d1 = target.raw_id();
                    component->detectors = {d1};
                } else if (num_translated_relative_detectors == 2) {
                    component->detectors = {component->detectors.d1, target.raw_id()};
                } else {
                    // We mark errors which have more than 3 detectors as a special boundary-to-boundary edge.
                    component->detectors = {SIZE_MAX, SIZE_MAX};
                }
            } else if (target.is_observable_id()) {
                component->observables.push_back(target.val());
            } else if (target.is_separator()) {
                // If the previous error in the decomposition had more than 3 components, we ignore it.
                if (!decomposed_err.components.back().detectors.is_valid_edge()) {
                    decomposed_err.components.pop_back();
                } else if (p > 0) {
                    handle_dem_error(p, component->detectors, component->observables);
                }
                decomposed_err.components.push_back({});
                component = &decomposed_err.components.back();
                num_translated_relative_detectors = 0;
            }
        }
        // If the final error in the decomposition had more than 3 components, we ignore it.
        if (!decomposed_err.components.back().detectors.is_valid_edge()) {
            decomposed_err.components.pop_back();
        } else if (p > 0) {
            handle_dem_error(p, component->detectors, component->observables);
        }

        // Add error to conditional groups.
        add_decomposed_error_to_conditional_groups(decomposed_err, conditional_groups);
    });
}

MatchingGraph detector_error_model_to_matching_graph(
    const stim::DetectorErrorModel &detector_error_model, pm::weight_int num_distinct_weights);

}  // namespace pm

#endif  // PYMATCHING2_IO_H
