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

#ifndef PYMATCHING2_USER_GRAPH_H
#define PYMATCHING2_USER_GRAPH_H

#include <cmath>
#include <list>
#include <set>
#include <stdexcept>
#include <vector>

#include "pymatching/sparse_blossom/driver/implied_weights.h"
#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/ints.h"
#include "pymatching/sparse_blossom/matcher/mwpm.h"
#include "pymatching/sparse_blossom/search/search_graph.h"
#include "stim.h"

namespace pm {

const pm::weight_int NUM_DISTINCT_WEIGHTS = 1 << (sizeof(pm::weight_int) * 8 - 8);

struct UserEdge {
    size_t node1;
    size_t node2;
    std::vector<size_t> observable_indices;  /// The indices of the observables crossed along this edge
    double weight;                           /// The weight of the edge to this neighboring node
    double error_probability;                /// The error probability associated with this node
    std::vector<ImpliedWeightUnconverted> implied_weights_for_other_edges;
};

struct UserNeighbor {
    std::list<UserEdge>::iterator edge_it;
    uint8_t pos{};  // The position of the neighboring node in the edge (either 0 if node1, or 1 if node2)
};

class UserNode {
   public:
    UserNode();
    size_t index_of_neighbor(size_t node) const;
    std::vector<UserNeighbor> neighbors;  /// The node's neighbors.
    bool is_boundary;
};

const pm::weight_int MAX_USER_EDGE_WEIGHT = NUM_DISTINCT_WEIGHTS - 1;

enum MERGE_STRATEGY : uint8_t { DISALLOW, INDEPENDENT, SMALLEST_WEIGHT, KEEP_ORIGINAL, REPLACE };

class UserGraph {
   public:
    std::vector<UserNode> nodes;
    std::list<UserEdge> edges;
    std::set<size_t> boundary_nodes;
    bool loaded_from_dem_without_correlations = false;

    UserGraph();
    explicit UserGraph(size_t num_nodes);
    UserGraph(size_t num_nodes, size_t num_observables);
    void merge_edge_or_boundary_edge(
        size_t node,
        size_t neighbor_index,
        const std::vector<size_t>& parallel_observables,
        double parallel_weight,
        double parallel_error_probability,
        MERGE_STRATEGY merge_strategy = DISALLOW);
    void add_or_merge_edge(
        size_t node1,
        size_t node2,
        const std::vector<size_t>& observables,
        double weight,
        double error_probability,
        MERGE_STRATEGY merge_strategy = DISALLOW);
    void add_or_merge_boundary_edge(
        size_t node,
        const std::vector<size_t>& observables,
        double weight,
        double error_probability,
        MERGE_STRATEGY merge_strategy = DISALLOW);
    bool has_edge(size_t node1, size_t node2);
    bool has_boundary_edge(size_t node);
    bool get_edge_or_boundary_edge_weight(size_t node1, size_t node2, double& weight_out);
    void set_boundary(const std::set<size_t>& boundary);
    std::set<size_t> get_boundary();
    size_t get_num_observables();
    void set_min_num_observables(size_t num_observables);
    size_t get_num_nodes();
    size_t get_num_detectors();
    size_t get_num_edges();
    bool is_boundary_node(size_t node_id);
    void add_noise(uint8_t* error_arr, uint8_t* syndrome_arr) const;
    bool all_edges_have_error_probabilities();
    double max_abs_weight();
    double get_edge_weight_normalising_constant(size_t max_num_distinct_weights);
    template <typename EdgeCallable, typename BoundaryEdgeCallable>
    double iter_discretized_edges(
        pm::weight_int num_distinct_weights,
        const EdgeCallable& edge_func,
        const BoundaryEdgeCallable& boundary_edge_func);
    template <typename EdgeCallable, typename BoundaryEdgeCallable>
    double to_matching_or_search_graph_helper(
        pm::weight_int num_distinct_weights,
        const EdgeCallable& edge_func,
        const BoundaryEdgeCallable& boundary_edge_func);
    pm::MatchingGraph to_matching_graph(pm::weight_int num_distinct_weights);
    pm::SearchGraph to_search_graph(pm::weight_int num_distinct_weights);
    pm::Mwpm to_mwpm(pm::weight_int num_distinct_weights, bool ensure_search_graph_included);
    void update_mwpm();
    Mwpm& get_mwpm();
    Mwpm& get_mwpm_with_search_graph();
    void handle_dem_instruction(double p, const std::vector<size_t>& detectors, const std::vector<size_t>& observables);
    void handle_dem_instruction_include_correlations(
        double p, const std::vector<size_t>& detectors, const std::vector<size_t>& observables);
    void get_nodes_on_shortest_path_from_source(size_t src, size_t dst, std::vector<size_t>& out_nodes);
    void populate_implied_edge_weights(
        std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>>& joint_probabilites);

   private:
    pm::Mwpm _mwpm;
    size_t _num_observables;
    bool _mwpm_needs_updating;
    bool _all_edges_have_error_probabilities;
};

double to_weight_for_correlations(double probability);

template <typename EdgeCallable, typename BoundaryEdgeCallable>
inline double UserGraph::iter_discretized_edges(
    pm::weight_int num_distinct_weights,
    const EdgeCallable& edge_func,
    const BoundaryEdgeCallable& boundary_edge_func) {
    double normalising_constant = get_edge_weight_normalising_constant(num_distinct_weights);

    for (auto& e : edges) {
        pm::signed_weight_int w = (pm::signed_weight_int)round(e.weight * normalising_constant);
        // Extremely important!
        // If all edge weights are even integers, then all collision events occur at integer times.
        w *= 2;
        bool node1_boundary = is_boundary_node(e.node1);
        bool node2_boundary = is_boundary_node(e.node2);
        if (node2_boundary && !node1_boundary) {
            boundary_edge_func(e.node1, w, e.observable_indices, e.implied_weights_for_other_edges);
        } else if (node1_boundary && !node2_boundary) {
            boundary_edge_func(e.node2, w, e.observable_indices, e.implied_weights_for_other_edges);
        } else if (!node1_boundary) {
            edge_func(e.node1, e.node2, w, e.observable_indices, e.implied_weights_for_other_edges);
        }
    }
    return normalising_constant * 2;
}

template <typename EdgeCallable, typename BoundaryEdgeCallable>
inline double UserGraph::to_matching_or_search_graph_helper(
    pm::weight_int num_distinct_weights,
    const EdgeCallable& edge_func,
    const BoundaryEdgeCallable& boundary_edge_func) {
    // Use vectors to store boundary edges initially before adding them to the graph, so
    // that parallel boundary edges with negative edge weights can be handled correctly
    std::vector<bool> has_boundary_edge(nodes.size(), false);
    std::vector<pm::signed_weight_int> boundary_edge_weights(nodes.size());
    std::vector<std::vector<size_t>> boundary_edge_observables(nodes.size());
    std::vector<std::vector<ImpliedWeightUnconverted>> boundary_edge_implied_weights_unconverted(nodes.size());

    double normalising_constant = iter_discretized_edges(
        num_distinct_weights,
        edge_func,
        [&](size_t u,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_for_other_edges) {
            // For parallel boundary edges, keep the boundary edge with the smaller weight
            if (!has_boundary_edge[u] || boundary_edge_weights[u] > weight) {
                boundary_edge_weights[u] = weight;
                boundary_edge_observables[u] = observables;
                has_boundary_edge[u] = true;
                boundary_edge_implied_weights_unconverted[u] = implied_weights_for_other_edges;
            }
        });

    // Now add boundary edges to the graph
    for (size_t i = 0; i < has_boundary_edge.size(); i++) {
        if (has_boundary_edge[i])
            boundary_edge_func(
                i,
                boundary_edge_weights[i],
                boundary_edge_observables[i],
                boundary_edge_implied_weights_unconverted[i]);
    }
    return normalising_constant;
}

UserGraph detector_error_model_to_user_graph(
    const stim::DetectorErrorModel& detector_error_model,
    bool enable_correlations,
    pm::weight_int num_distinct_weights);

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
    const stim::DetectorErrorModel& detector_error_model, const Handler& handle_dem_error) {
    detector_error_model.iter_flatten_error_instructions([&](const stim::DemInstruction& instruction) {
        std::vector<size_t> dets;
        std::vector<size_t> observables;
        double p = instruction.arg_data[0];
        for (auto& target : instruction.target_data) {
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

struct DecomposedDemError {
    /// The probability of this error occurring.
    double probability;
    /// Effects of the error.
    std::vector<UserEdge> components;

    bool operator==(const DecomposedDemError& other) const;
    bool operator!=(const DecomposedDemError& other) const;
};

void add_decomposed_error_to_joint_probabilities(
    DecomposedDemError& error,
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>>& joint_probabilites);

template <typename Handler>
void iter_dem_instructions_include_correlations(
    const stim::DetectorErrorModel& detector_error_model,
    const Handler& handle_dem_error,
    std::map<std::pair<size_t, size_t>, std::map<std::pair<size_t, size_t>, double>>& joint_probabilites,
    bool include_decomposed_error_components_in_edge_weights = true) {
    detector_error_model.iter_flatten_error_instructions([&](const stim::DemInstruction& instruction) {
        double p = instruction.arg_data[0];
        pm::DecomposedDemError decomposed_err;
        decomposed_err.probability = p;
        if (p > 0.5) {
            throw ::std::invalid_argument(
                "Errors with probability greater than 0.5 are not supported with correlations enabled");
        }
        if (p == 0) {
            // Ignore errors with no error probability.
            return;
        }
        decomposed_err.components = {};
        decomposed_err.components.push_back({});
        UserEdge* component = &decomposed_err.components.back();
        // Mark component as empty to begin with.
        component->node1 = SIZE_MAX;
        component->node2 = SIZE_MAX;
        size_t num_component_detectors = 0;
        bool instruction_contains_separator = false;
        for (auto& target : instruction.target_data) {
            // Decompose error
            if (target.is_relative_detector_id()) {
                num_component_detectors++;
                if (num_component_detectors == 1) {
                    const size_t& d1 = target.raw_id();
                    component->node1 = d1;
                } else if (num_component_detectors == 2) {
                    // Maintain invariant that node1 <= node2.
                    if (component->node1 <= target.raw_id()) {
                        component->node2 = target.raw_id();
                    } else {
                        component->node2 = component->node1;
                        component->node1 = target.raw_id();
                    }
                } else {
                    // Undecomposed hyperedges are not supported
                    throw std::invalid_argument(
                        "Encountered an undecomposed error instruction with 3 or mode detectors. "
                        "This is not supported when using `enable_correlations=True`. "
                        "Did you forget to set `decompose_errors=True` when "
                        "converting the stim circuit to a detector error model?");
                }
            } else if (target.is_observable_id()) {
                component->observable_indices.push_back(target.val());
            } else if (target.is_separator()) {
                instruction_contains_separator = true;
                // We cannot have num_component_detectors > 2 at this point, or we would have already thrown an
                // exception
                if (num_component_detectors == 0) {
                    throw std::invalid_argument(
                        "Encountered a decomposed error instruction with an undetectable component (0 detectors). "
                        "This is not supported.");
                }
                // The previous error in the decomposition must have 1 or 2 detectors
                decomposed_err.components.push_back({});
                component = &decomposed_err.components.back();
                component->node1 = SIZE_MAX;
                component->node2 = SIZE_MAX;
                num_component_detectors = 0;
            }
        }

        if (num_component_detectors == 0) {
            if (instruction_contains_separator) {
                throw std::invalid_argument(
                    "Encountered a decomposed error instruction with an undetectable component (0 detectors). "
                    "This is not supported.");
            } else {
                // Ignore errors that are undetectable, provided they are not a component of a decomposed error
                return;
            }
        }

        // If include_decomposed_error_components_in_edge_weights is False, then only add the edge into 
        // the graph if it is not a component in a decomposed error with more than one component
        if (include_decomposed_error_components_in_edge_weights || decomposed_err.components.size() == 1) {
            for (pm::UserEdge& component : decomposed_err.components) {
                if (component.node2 == SIZE_MAX) {
                    handle_dem_error(p, {component.node1}, component.observable_indices);
                } else {
                    handle_dem_error(p, {component.node1, component.node2}, component.observable_indices);
                }
            }
        }

        add_decomposed_error_to_joint_probabilities(decomposed_err, joint_probabilites);
    });
}

}  // namespace pm

#endif  // PYMATCHING2_USER_GRAPH_H
