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

#include <algorithm>
#include <utility>

#include "pymatching/sparse_blossom/flooder/graph.h"

double pm::merge_weights(double a, double b) {
    auto sgn = std::copysign(1, a) * std::copysign(1, b);
    auto signed_min = sgn * std::min(std::abs(a), std::abs(b));
    return signed_min + std::log(1 + std::exp(-std::abs(a + b))) - std::log(1 + std::exp(-std::abs(a - b)));
}

void pm::IntermediateWeightedGraph::add_or_merge_edge(
    size_t u, size_t v, double weight, const std::vector<size_t>& observables) {
    size_t larger_node = std::max(u, v);
    if (larger_node + 1 > nodes.size()) {
        throw std::invalid_argument(
            "Node " + std::to_string(larger_node) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }
    auto it = std::find_if(nodes[u].begin(), nodes[u].end(), [this, v](Neighbor& neighbor) {
        return neighbor.node == &nodes[v];
    });
    if (it == nodes[u].end()) {
        nodes[u].push_back({&nodes[v], weight, observables});
        nodes[v].push_back({&nodes[u], weight, observables});
    } else {
        double new_weight = merge_weights((*it).weight, weight);
        (*it).weight = new_weight;
        auto it2 = std::find_if((*it).node->begin(), (*it).node->end(), [this, u](Neighbor& neighbor) {
            return neighbor.node == &nodes[u];
        });
        if (it2 != (*it).node->end())
            (*it2).weight = new_weight;
    }
}

void pm::IntermediateWeightedGraph::add_or_merge_boundary_edge(
    size_t u, double weight, const std::vector<size_t>& observables) {
    if (u > nodes.size() - 1) {
        throw std::invalid_argument(
            "Node " + std::to_string(u) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }
    auto it = std::find_if(nodes[u].begin(), nodes[u].end(), [this](Neighbor& neighbor) {
        return neighbor.node == nullptr;
    });
    if (it == nodes[u].end()) {
        nodes[u].push_back({nullptr, weight, observables});
    } else {
        (*it).weight = merge_weights(weight, (*it).weight);
    }
}

void pm::IntermediateWeightedGraph::handle_dem_instruction(
    double p, const std::vector<size_t>& detectors, const std::vector<size_t>& observables) {
    if (detectors.size() == 2) {
        add_or_merge_edge(detectors[0], detectors[1], std::log((1 - p) / p), observables);
    } else if (detectors.size() == 1) {
        add_or_merge_boundary_edge(detectors[0], std::log((1 - p) / p), observables);
    }
}

double pm::IntermediateWeightedGraph::max_abs_weight() {
    double max_abs_weight = 0;
    for (auto& node : nodes) {
        for (auto& neighbor : node) {
            if (std::abs(neighbor.weight) > max_abs_weight)
                max_abs_weight = std::abs(neighbor.weight);
        }
    }
    return max_abs_weight;
}

pm::MatchingGraph pm::IntermediateWeightedGraph::to_matching_graph(pm::weight_int num_distinct_weights) {
    pm::MatchingGraph matching_graph(nodes.size(), num_observables);
    double normalising_constant = iter_discretized_edges(
        num_distinct_weights,
        [&](size_t u, size_t v, pm::signed_weight_int weight, const std::vector<size_t>& observables) {
            matching_graph.add_edge(u, v, weight, observables);
        },
        [&](size_t u, pm::signed_weight_int weight, const std::vector<size_t>& observables) {
            matching_graph.add_boundary_edge(u, weight, observables);
        });
    matching_graph.normalising_constant = normalising_constant;
    return matching_graph;
}

pm::SearchGraph pm::IntermediateWeightedGraph::to_search_graph(pm::weight_int num_distinct_weights) {
    /// Identical to pm::IntermediateWeightedGraph::to_matching_graph but for constructing a pm::SearchGraph
    pm::SearchGraph search_graph(nodes.size());
    iter_discretized_edges(
        num_distinct_weights,
        [&](size_t u, size_t v, pm::signed_weight_int weight, const std::vector<size_t>& observables) {
            search_graph.add_edge(u, v, weight, observables);
        },
        [&](size_t u, pm::signed_weight_int weight, const std::vector<size_t>& observables) {
            search_graph.add_boundary_edge(u, weight, observables);
        });
    return search_graph;
}

pm::Mwpm pm::IntermediateWeightedGraph::to_mwpm(
    pm::weight_int num_distinct_weights, bool ensure_search_flooder_included) {
    if (num_observables > sizeof(pm::obs_int) * 8 || ensure_search_flooder_included) {
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

pm::IntermediateWeightedGraph pm::detector_error_model_to_weighted_graph(
    const stim::DetectorErrorModel& detector_error_model) {
    pm::IntermediateWeightedGraph weighted_graph(
        detector_error_model.count_detectors(), detector_error_model.count_observables());
    pm::iter_detector_error_model_edges(
        detector_error_model, [&](double p, const std::vector<size_t>& detectors, std::vector<size_t>& observables) {
            weighted_graph.handle_dem_instruction(p, detectors, observables);
        });
    return weighted_graph;
}

pm::MatchingGraph pm::detector_error_model_to_matching_graph(
    const stim::DetectorErrorModel& detector_error_model, pm::weight_int num_distinct_weights) {
    auto weighted_graph = pm::detector_error_model_to_weighted_graph(detector_error_model);
    return weighted_graph.to_matching_graph(num_distinct_weights);
}
bool pm::Neighbor::operator==(const pm::Neighbor& rhs) const {
    return node == rhs.node && weight == rhs.weight && observables == rhs.observables;
}
bool pm::Neighbor::operator!=(const pm::Neighbor& rhs) const {
    return !(rhs == *this);
}
