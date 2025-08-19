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

#include "search_graph.h"

#include <cassert>
#include <cmath>

pm::SearchGraph::SearchGraph() : num_nodes(0) {
}

pm::SearchGraph::SearchGraph(size_t num_nodes) : num_nodes(num_nodes) {
    nodes.resize(num_nodes);
}

pm::SearchGraph::SearchGraph(pm::SearchGraph&& graph) noexcept
    : nodes(std::move(graph.nodes)),
      num_nodes(graph.num_nodes),
      negative_weight_edges(std::move(graph.negative_weight_edges)) {
}

void pm::SearchGraph::add_edge(
    size_t u,
    size_t v,
    signed_weight_int weight,
    const std::vector<size_t>& observables,
    const std::vector<ImpliedWeightUnconverted>& implied_weights_for_other_edges) {
    size_t larger_node = std::max(u, v);
    if (larger_node + 1 > nodes.size()) {
        throw std::invalid_argument(
            "Node " + std::to_string(larger_node) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }

    // Self-loops are ignored, since they are never included in the shortest-path (the search graph has only
    // positive weights).
    if (u == v)
        return;

    uint8_t weight_sign = 0;
    if (weight < 0) {
        negative_weight_edges.push_back({u, v});
        weight_sign = pm::WEIGHT_SIGN;
    }

    nodes[u].neighbors.push_back(&(nodes[v]));
    nodes[u].neighbor_weights.push_back(std::abs(weight));
    nodes[u].neighbor_observable_indices.push_back(observables);
    nodes[u].neighbor_markers.push_back(weight_sign);
    nodes[u].neighbor_implied_weights.push_back({});
    edges_to_implied_weights_unconverted[u].emplace_back(implied_weights_for_other_edges);

    nodes[v].neighbors.push_back(&(nodes[u]));
    nodes[v].neighbor_weights.push_back(std::abs(weight));
    nodes[v].neighbor_observable_indices.push_back(observables);
    nodes[v].neighbor_markers.push_back(weight_sign);
    nodes[v].neighbor_implied_weights.push_back({});
    edges_to_implied_weights_unconverted[v].emplace_back(implied_weights_for_other_edges);
}

void pm::SearchGraph::add_boundary_edge(
    size_t u,
    signed_weight_int weight,
    const std::vector<size_t>& observables,
    const std::vector<ImpliedWeightUnconverted>& implied_weights_for_other_edges) {
    if (u >= nodes.size()) {
        throw std::invalid_argument(
            "Node " + std::to_string(u) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }

    uint8_t weight_sign = 0;
    if (weight < 0) {
        negative_weight_edges.push_back({u, SIZE_MAX});
        weight_sign = pm::WEIGHT_SIGN;
    }

    nodes[u].neighbors.insert(nodes[u].neighbors.begin(), 1, nullptr);
    nodes[u].neighbor_weights.insert(nodes[u].neighbor_weights.begin(), 1, std::abs(weight));
    nodes[u].neighbor_observable_indices.insert(nodes[u].neighbor_observable_indices.begin(), 1, observables);
    nodes[u].neighbor_markers.insert(nodes[u].neighbor_markers.begin(), 1, weight_sign);
    nodes[u].neighbor_implied_weights.insert(nodes[u].neighbor_implied_weights.begin(), 1, {});
    edges_to_implied_weights_unconverted[u].insert(
        edges_to_implied_weights_unconverted[u].begin(), 1, implied_weights_for_other_edges);
}

// Reweight assuming an error has occurred on a single edge u, v. When v == -1, assumes an edge from
// u to the boundary.
void pm::SearchGraph::reweight_for_edge(const int64_t& u, const int64_t& v) {
    size_t z = nodes[u].index_of_neighbor(v == -1 ? nullptr : &nodes[v]);
    reweight(nodes[u].neighbor_implied_weights[z]);
}

void pm::SearchGraph::reweight_for_edges(const std::vector<int64_t>& edges) {
    for (size_t i = 0; i < edges.size() >> 1; ++i) {
        int64_t u = edges[2 * i];
        int64_t v = edges[2 * i + 1];
        reweight_for_edge(u, v);
    }
}

void pm::SearchGraph::undo_reweights() {
    // We iterate backward over the previous weights, since some edges
    // may have been reweighted more than once.
    for (auto it = previous_weights.rbegin(); it != previous_weights.rend(); ++it) {
        pm::PreviousWeight& prev = *it;
        *prev.ptr = prev.val;
    }
    previous_weights.clear();
}

namespace {

pm::ImpliedWeight convert_rule(
    std::vector<pm::SearchDetectorNode>& nodes,
    const pm::ImpliedWeightUnconverted& rule,
    const double normalising_constant) {
    const size_t& i = rule.node1;
    const size_t& j = rule.node2;
    pm::weight_int* weight_pointer_i =
        &nodes[i].neighbor_weights[nodes[i].index_of_neighbor(j == SIZE_MAX ? nullptr : &nodes[j])];
    pm::weight_int* weight_pointer_j =
        j == SIZE_MAX ? nullptr : &nodes[j].neighbor_weights[nodes[j].index_of_neighbor(&nodes[i])];

    double rescaled_normalising_constant = normalising_constant / 2;
    pm::signed_weight_int w = (pm::signed_weight_int)round(rule.implied_weight * rescaled_normalising_constant);
    // Extremely important!
    // If all edge weights are even integers, then all collision events occur at integer times.
    w *= 2;
    return pm::ImpliedWeight{weight_pointer_i, weight_pointer_j, static_cast<pm::weight_int>(std::abs(w))};
}

}  // namespace

void pm::SearchGraph::convert_implied_weights(const double normalising_constant) {
    for (size_t u = 0; u < nodes.size(); u++) {
        const std::vector<std::vector<ImpliedWeightUnconverted>>& rules_for_node =
            edges_to_implied_weights_unconverted[u];
        for (size_t v = 0; v < nodes[u].neighbors.size(); v++) {
            for (const auto& rule : rules_for_node[v]) {
                ImpliedWeight converted = convert_rule(nodes, rule, normalising_constant);
                nodes[u].neighbor_implied_weights[v].push_back(converted);
            }
        }
    }
}
