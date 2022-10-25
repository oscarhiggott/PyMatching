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
#include "pymatching/sparse_blossom/matcher/mwpm.h"
#include "pymatching/sparse_blossom/search/search_graph.h"
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

struct Neighbor {
    std::vector<Neighbor> *node;
    bool operator==(const Neighbor &rhs) const;
    bool operator!=(const Neighbor &rhs) const;
    double weight;
    std::vector<size_t> observables;
};

class IntermediateWeightedGraph {
   public:
    std::vector<std::vector<Neighbor>> nodes;
    size_t num_nodes;
    size_t num_observables;

    explicit IntermediateWeightedGraph(size_t num_nodes, size_t num_observables)
        : num_nodes(num_nodes), num_observables(num_observables) {
        nodes.resize(num_nodes);
    };

    void add_or_merge_edge(size_t u, size_t v, double weight, const std::vector<size_t> &observables);

    void add_or_merge_boundary_edge(size_t u, double weight, const std::vector<size_t> &observables);

    void handle_dem_instruction(double p, const std::vector<size_t> &detectors, const std::vector<size_t> &observables);

    template <typename EdgeCallable, typename BoundaryEdgeCallable>
    double iter_discretized_edges(
        pm::weight_int num_distinct_weights,
        const EdgeCallable &edge_func,
        const BoundaryEdgeCallable &boundary_edge_func);

    pm::MatchingGraph to_matching_graph(pm::weight_int num_distinct_weights);

    pm::SearchGraph to_search_graph(pm::weight_int num_distinct_weights);

    pm::Mwpm to_mwpm(pm::weight_int num_distinct_weights, bool ensure_search_flooder_included = false);

    double max_abs_weight();
};

template <typename EdgeCallable, typename BoundaryEdgeCallable>
inline double IntermediateWeightedGraph::iter_discretized_edges(
    pm::weight_int num_distinct_weights,
    const EdgeCallable &edge_func,
    const BoundaryEdgeCallable &boundary_edge_func) {
    double max_weight = max_abs_weight();
    pm::MatchingGraph matching_graph(nodes.size(), num_observables);
    pm::weight_int max_half_edge_weight = num_distinct_weights - 1;
    double normalising_constant = (double)max_half_edge_weight / max_weight;
    for (auto &node : nodes) {
        for (auto &neighbor : node) {
            auto i = &node - &nodes[0];
            pm::signed_weight_int w = (pm::signed_weight_int)round(neighbor.weight * normalising_constant);

            // Extremely important!
            // If all edge weights are even integers, then all collision events occur at integer times.
            w *= 2;

            if (!neighbor.node) {
                boundary_edge_func(i, w, neighbor.observables);
            } else {
                auto j = neighbor.node - &nodes[0];
                if (j > i)
                    edge_func(i, j, w, neighbor.observables);
            }
        }
    }
    return normalising_constant * 2;
}

IntermediateWeightedGraph detector_error_model_to_weighted_graph(const stim::DetectorErrorModel &detector_error_model);

MatchingGraph detector_error_model_to_matching_graph(
    const stim::DetectorErrorModel &detector_error_model, pm::weight_int num_distinct_weights);

}  // namespace pm

#endif  // PYMATCHING2_IO_H
