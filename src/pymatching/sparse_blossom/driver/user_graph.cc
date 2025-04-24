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

#include <cstdint>
#include <pymatching/sparse_blossom/flooder/detector_node.h>

#include "pymatching/rand/rand_gen.h"
#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/matcher/mwpm.h"
#include "pymatching/sparse_blossom/search/search_graph.h"

bool pm::UserGraph::contains(size_t node1, size_t node2) {
    if (edge_map.count({node1, node2})) {
        return true;
    }
    return false;
}

std::optional<pm::DetectorEdgeData> pm::UserGraph::get_edge_data(size_t node1, size_t node2) {
    DetectorEdgeId search_node{node1, node2};
    if (contains(node1, node2)) {
        return edge_map[search_node];
    }
    return std::nullopt;
}

bool is_valid_probability(double p) {
    return (p >= 0 && p <= 1);
}

void pm::UserGraph::merge_edge_or_boundary_edge(
    size_t node1,
    size_t node2,
    const std::vector<size_t>& new_observables,
    double cur_edge_weight,
    double parallel_edge_weight,
    double cur_edge_probability,
    double parallel_edge_probability,
    pm::MERGE_STRATEGY merge_strategy) {
    DetectorEdgeId existing_edge{node1, node2};
    if (merge_strategy == DISALLOW) {
        throw std::invalid_argument(
            "Edge (" + std::to_string(existing_edge.d1) + ", " + std::to_string(existing_edge.d2) +
            ") already exists in the graph. "
            "Parallel edges not permitted with the provided `disallow` `merge_strategy`. Please provide a "
            "different `merge_strategy`.");
    } else if (
        merge_strategy == KEEP_ORIGINAL ||
        (merge_strategy == SMALLEST_WEIGHT && parallel_edge_weight >= cur_edge_weight)) {
        return;
    } else {
        double new_weight, new_error_probability;
        bool use_new_observables;
        if (merge_strategy == REPLACE || merge_strategy == SMALLEST_WEIGHT) {
            new_weight = parallel_edge_weight;
            new_error_probability = parallel_edge_probability;
            use_new_observables = true;
        } else if (merge_strategy == INDEPENDENT) {
            new_weight = pm::merge_weights(cur_edge_weight, parallel_edge_weight);
            new_error_probability = -1;
            if (is_valid_probability(cur_edge_probability) && is_valid_probability(parallel_edge_probability))
                new_error_probability = parallel_edge_probability * (1 - cur_edge_probability) +
                                        cur_edge_probability * (1 - parallel_edge_probability);
            // We do not need to update the observables. If they do not match up, then the code has distance 2.
            use_new_observables = false;
        } else {
            throw std::invalid_argument("Merge strategy not recognised.");
        }
        // Update the existing edge weight and probability in the adjacency list of `node`
        edge_map[existing_edge].error_probability = new_error_probability;
        edge_map[existing_edge].weight = new_weight;
        if (use_new_observables) {
            edge_map[existing_edge].observables = new_observables;
        }

        _mwpm_needs_updating = true;
        if (new_error_probability < 0 || new_error_probability > 1)
            _all_edges_have_error_probabilities = false;
    }
}

void pm::UserGraph::add_or_merge_edge(
    size_t node1,
    size_t node2,
    const std::vector<size_t>& observables,
    double new_weight,
    double new_error_probability,
    MERGE_STRATEGY merge_strategy) {
    auto max_id = std::max(node1, node2);
    if (max_id + 1 > num_nodes) {
        num_nodes = max_id + 1;
    }

    bool edge_exists = contains(node1, node2);

    if (!edge_exists) {
        pm::DetectorEdgeData edge = {
            .detectors = {node1, node2},
            .observables = observables,
            .weight = new_weight,
            .error_probability = new_error_probability,
            .correlated_proabilities_sum = 0,
            .implied_weights_for_other_edges = {}};
        edge_map[{node1, node2}] = edge;
        for (auto& obs : observables) {
            if (obs + 1 > _num_observables)
                _num_observables = obs + 1;
        }
        _mwpm_needs_updating = true;
        if (new_error_probability < 0 || new_error_probability > 1)
            _all_edges_have_error_probabilities = false;
    } else {
        DetectorEdgeData& edge = edge_map[{node1, node2}];
        merge_edge_or_boundary_edge(
            node1,
            node2,
            observables,
            edge.weight,
            new_weight,
            edge.error_probability,
            new_error_probability,
            merge_strategy);
    }
}

void pm::UserGraph::add_or_merge_boundary_edge(
    size_t node,
    const std::vector<size_t>& observables,
    double new_weight,
    double new_error_probability,
    MERGE_STRATEGY merge_strategy) {
    if (node + 1 > num_nodes) {
        num_nodes = node + 1;
    }

    bool edge_exists = contains(node, SIZE_MAX);

    if (!edge_exists) {
        pm::DetectorEdgeData edge = {{node, SIZE_MAX}, observables, new_weight, new_error_probability};
        edge_map[{node, SIZE_MAX}] = {
            .detectors = {node, SIZE_MAX},
            .observables = observables,
            .weight = new_weight,
            .error_probability = new_error_probability,
            .correlated_proabilities_sum = 0,
            .implied_weights_for_other_edges = {}};

        for (auto& obs : observables) {
            if (obs + 1 > _num_observables)
                _num_observables = obs + 1;
        }
        _mwpm_needs_updating = true;
        if (new_error_probability < 0 || new_error_probability > 1)
            _all_edges_have_error_probabilities = false;
    } else {
        DetectorEdgeData& edge = edge_map[{node, SIZE_MAX}];
        merge_edge_or_boundary_edge(
            node,
            SIZE_MAX,
            observables,
            edge.weight,
            new_weight,
            edge.error_probability,
            new_error_probability,
            merge_strategy);
    }
}

pm::UserGraph::UserGraph()
    : _num_observables(0), _mwpm_needs_updating(true), _all_edges_have_error_probabilities(true) {
}

pm::UserGraph::UserGraph(size_t num_nodes)
    : _num_observables(0), _mwpm_needs_updating(true), _all_edges_have_error_probabilities(true), num_nodes(num_nodes) {
}

pm::UserGraph::UserGraph(size_t num_nodes, size_t num_observables)
    : _num_observables(num_observables),
      _mwpm_needs_updating(true),
      _all_edges_have_error_probabilities(true),
      num_nodes(num_nodes) {
}

void pm::UserGraph::set_boundary(const std::set<size_t>& boundary) {
    boundary_nodes = std::move(boundary);
    for (auto& n : boundary_nodes) {
        if (n >= num_nodes)
            num_nodes = n;
    }
    _mwpm_needs_updating = true;
}

std::set<size_t> pm::UserGraph::get_boundary() {
    return boundary_nodes;
}

size_t pm::UserGraph::get_num_observables() {
    return _num_observables;
}

size_t pm::UserGraph::get_num_nodes() {
    return num_nodes;
}

size_t pm::UserGraph::get_num_detectors() {
    return get_num_nodes() - boundary_nodes.size();
}

bool pm::UserGraph::is_boundary_node(size_t node) {
    return (node == SIZE_MAX) || boundary_nodes.contains(node);
}

void pm::UserGraph::update_mwpm() {
    _mwpm = to_mwpm(pm::NUM_DISTINCT_WEIGHTS, false);
    _mwpm_needs_updating = false;
}

pm::Mwpm& pm::UserGraph::get_mwpm() {
    if (_mwpm_needs_updating)
        update_mwpm();
    return _mwpm;
}

void pm::UserGraph::add_noise(uint8_t* error_arr, uint8_t* syndrome_arr) const {
    if (!_all_edges_have_error_probabilities)
        return;

    for (auto& edgeIDandData : edge_map) {
        auto p = edgeIDandData.second.error_probability;
        if (rand_float(0.0, 1.0) < p) {
            // Flip the observables
            for (auto& obs : edgeIDandData.second.observables) {
                *(error_arr + obs) ^= 1;
            }
            // Flip the syndrome bits
            *(syndrome_arr + edgeIDandData.first.d1) ^= 1;
            if (edgeIDandData.first.d2 != SIZE_MAX)
                *(syndrome_arr + edgeIDandData.first.d2) ^= 1;
        }
    }

    for (auto& b : boundary_nodes)
        *(syndrome_arr + b) = 0;
}

size_t pm::UserGraph::get_num_edges() {
    return edge_map.size();
}

bool pm::UserGraph::all_edges_have_error_probabilities() {
    return _all_edges_have_error_probabilities;
}

double pm::UserGraph::max_abs_weight() {
    double max_abs_weight = 0;
    for (auto& edgeIdAndData : edge_map) {
        if (std::abs(edgeIdAndData.second.weight) > max_abs_weight) {
            max_abs_weight = std::abs(edgeIdAndData.second.weight);
        }
    }
    return max_abs_weight;
}

pm::MatchingGraph pm::UserGraph::to_matching_graph(pm::weight_int num_distinct_weights) {
    pm::MatchingGraph matching_graph(get_num_nodes(), _num_observables);

    size_t num_edges = 0;
    double normalising_constant = to_matching_or_search_graph_helper(
        num_distinct_weights,
        [&](size_t u,
            size_t v,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_unconverted) {
            num_edges++;
            matching_graph.add_edge(u, v, weight, observables, implied_weights_unconverted);
        },
        [&](size_t u,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_unconverted) {
            num_edges++;
            matching_graph.add_boundary_edge(u, weight, observables, implied_weights_unconverted);
        });

    matching_graph.normalising_constant = normalising_constant;
    if (boundary_nodes.size() > 0) {
        matching_graph.is_user_graph_boundary_node.clear();
        matching_graph.is_user_graph_boundary_node.resize(get_num_nodes(), false);
        for (auto& i : boundary_nodes)
            matching_graph.is_user_graph_boundary_node[i] = true;
    }
    matching_graph.convert_implied_weights();
    return matching_graph;
}

pm::SearchGraph pm::UserGraph::to_search_graph(pm::weight_int num_distinct_weights) {
    /// Identical to to_matching_graph but for constructing a pm::SearchGraph
    pm::SearchGraph search_graph(get_num_nodes());

    to_matching_or_search_graph_helper(
        num_distinct_weights,
        [&](size_t u,
            size_t v,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_unconverted) {
            search_graph.add_edge(u, v, weight, observables);
        },
        [&](size_t u,
            pm::signed_weight_int weight,
            const std::vector<size_t>& observables,
            const std::vector<ImpliedWeightUnconverted>& implied_weights_unconverted) {
            search_graph.add_boundary_edge(u, weight, observables);
        });
    return search_graph;
}

pm::Mwpm pm::UserGraph::to_mwpm(pm::weight_int num_distinct_weights, bool ensure_search_graph_included) {
    if (_num_observables > sizeof(pm::obs_int) * 8 || ensure_search_graph_included) {
        auto mwpm = pm::Mwpm(
            pm::GraphFlooder(to_matching_graph(num_distinct_weights)),
            pm::SearchFlooder(to_search_graph(num_distinct_weights)));
        mwpm.flooder.sync_negative_weight_observables_and_detection_events();
        return mwpm;
    } else {
        auto mwpm = pm::Mwpm(pm::GraphFlooder(to_matching_graph(num_distinct_weights)));
        mwpm.flooder.sync_negative_weight_observables_and_detection_events();
        return mwpm;
    }
}

pm::Mwpm& pm::UserGraph::get_mwpm_with_search_graph() {
    if (!_mwpm_needs_updating && _mwpm.flooder.graph.nodes.size() == _mwpm.search_flooder.graph.nodes.size()) {
        return _mwpm;
    } else {
        _mwpm = to_mwpm(pm::NUM_DISTINCT_WEIGHTS, true);
        _mwpm_needs_updating = false;
        return _mwpm;
    }
}

void pm::UserGraph::handle_dem_instruction(
    double p, const DetectorEdgeId& edge, const std::vector<size_t>& observables) {
    if (edge.is_valid_edge()) {
        if (edge.is_boundary()) {
            add_or_merge_boundary_edge(edge.d1, observables, std::log((1 - p) / p), p, INDEPENDENT);
        } else {
            add_or_merge_edge(edge.d1, edge.d2, observables, std::log((1 - p) / p), p, INDEPENDENT);
        }
    }
}

void pm::UserGraph::get_nodes_on_shortest_path_from_source(size_t src, size_t dst, std::vector<size_t>& out_nodes) {
    auto& mwpm = get_mwpm_with_search_graph();
    bool src_is_boundary = is_boundary_node(src);
    bool dst_is_boundary = is_boundary_node(dst);
    if (src != SIZE_MAX && src >= get_num_nodes())
        throw std::invalid_argument("node " + std::to_string(src) + " is not in the graph");
    if (dst != SIZE_MAX && dst >= get_num_nodes())
        throw std::invalid_argument("node " + std::to_string(dst) + " is not in the graph");
    if (!src_is_boundary) {
        size_t dst_tmp = dst_is_boundary ? SIZE_MAX : dst;
        mwpm.search_flooder.iter_edges_on_shortest_path_from_source(src, dst_tmp, [&](const pm::SearchGraphEdge edge) {
            out_nodes.push_back(edge.detector_node - &mwpm.search_flooder.graph.nodes[0]);
        });
        if (!dst_is_boundary)
            out_nodes.push_back(dst);
    } else if (!dst_is_boundary) {
        std::vector<size_t> temp_out_nodes;
        get_nodes_on_shortest_path_from_source(dst, src, temp_out_nodes);
        for (size_t i = 0; i < temp_out_nodes.size(); i++) {
            out_nodes.push_back(temp_out_nodes[temp_out_nodes.size() - 1 - i]);
        }
    } else {
        throw std::invalid_argument("Both the source and destination vertices provided are boundary nodes");
    }
}

bool pm::UserGraph::has_edge(size_t node1, size_t node2) {
    if (node1 >= get_num_nodes())
        return false;
    return contains(node1, node2);
}

bool pm::UserGraph::has_boundary_edge(size_t node) {
    if (node >= get_num_nodes())
        return false;
    return contains(node, SIZE_MAX);
}

void pm::UserGraph::set_min_num_observables(size_t num_observables) {
    if (num_observables > _num_observables)
        _num_observables = num_observables;
}

double pm::UserGraph::get_edge_weight_normalising_constant(size_t max_num_distinct_weights) {
    double max_abs_weight = 0;
    bool all_integral_weight = true;
    for (auto& edgeIdAndData : edge_map) {
        DetectorEdgeData& e = edgeIdAndData.second;
        if (std::abs(e.weight) > max_abs_weight)
            max_abs_weight = std::abs(e.weight);

        if (round(e.weight) != e.weight)
            all_integral_weight = false;
    }

    if (max_abs_weight > pm::MAX_USER_EDGE_WEIGHT)
        throw std::invalid_argument(
            "maximum absolute edge weight of " + std::to_string(pm::MAX_USER_EDGE_WEIGHT) + " exceeded.");

    if (all_integral_weight) {
        return 1.0;
    } else {
        pm::weight_int max_half_edge_weight = max_num_distinct_weights - 1;
        return (double)max_half_edge_weight / max_abs_weight;
    }
}

pm::weight_int pm::convert_probability_to_weight(double p) {
    return std::log((1 - p) / p);
}

pm::UserGraph pm::detector_error_model_to_user_graph(const stim::DetectorErrorModel& detector_error_model) {
    pm::UserGraph user_graph(detector_error_model.count_detectors(), detector_error_model.count_observables());
    std::map<pm::DetectorEdgeId, std::map<pm::DetectorEdgeId, double>> conditional_groups;

    // Create initial user graph. Also calculate conditional_groups.
    pm::iter_dem_instructions_include_correlations(
        detector_error_model,
        [&](double p, const DetectorEdgeId& nodes, std::vector<size_t>& observables) {
            user_graph.handle_dem_instruction(p, nodes, observables);
        },
        conditional_groups);

    // Include reweighting information in the user graph.
    for (const auto& pf : conditional_groups) {
        auto causal_edge = pf.first;
        double summed_probabilities = 0;
        for (const auto& affected_edge_and_probability : pf.second) {
            summed_probabilities += affected_edge_and_probability.second;
        }
        user_graph.edge_map[causal_edge].correlated_proabilities_sum = summed_probabilities;
        for (const auto& affected_edge_and_probability : pf.second) {
            pm::DetectorEdgeId affected_edge = affected_edge_and_probability.first;
            if (affected_edge != causal_edge) {
                double implied_probability_for_other_edge =
                    std::min(0.9, affected_edge_and_probability.second / summed_probabilities);
                weight_int implied_weight_for_other_edge =
                    convert_probability_to_weight(implied_probability_for_other_edge);
                ImpliedWeightUnconverted iwu{affected_edge.d1, affected_edge.d2, implied_weight_for_other_edge};
                user_graph.edge_map[causal_edge].implied_weights_for_other_edges.emplace_back(iwu);
            }
        }
    }
    return user_graph;
}

void pm::UserGraph::set_num_nodes(const size_t n) {
    num_nodes = n;
}
