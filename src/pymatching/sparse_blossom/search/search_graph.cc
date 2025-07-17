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

void pm::SearchGraph::add_edge(size_t u, size_t v, signed_weight_int weight, const std::vector<size_t>& observables) {
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

    nodes[v].neighbors.push_back(&(nodes[u]));
    nodes[v].neighbor_weights.push_back(std::abs(weight));
    nodes[v].neighbor_observable_indices.push_back(observables);
    nodes[v].neighbor_markers.push_back(weight_sign);
}

void pm::SearchGraph::add_boundary_edge(size_t u, signed_weight_int weight, const std::vector<size_t>& observables) {
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
}

void pm::SearchGraph::reweight_for_edges(
    const std::vector<int64_t>& edges,
    ::pm::MatchingGraph& flooder_graph,
    const std::unordered_map<const ::pm::weight_int*, std::pair<int32_t, int32_t>>& flooder_graph_weight_ptr_map) {
    assert(previous_weights.empty());
    for (size_t i = 0; i < edges.size() >> 1; ++i) {
        int64_t u = edges[2 * i];
        int64_t v = edges[2 * i + 1];
        size_t z = flooder_graph.nodes[u].index_of_neighbor(v == -1 ? nullptr : &flooder_graph.nodes[v]);
        for (const auto& [flooder_edge0_ptr, flooder_edge1_ptr, implied_weight] :
             flooder_graph.nodes[u].neighbor_implied_weights[z]) {
            ::pm::weight_int* search_edge0_ptr = nullptr;
            ::pm::weight_int* search_edge1_ptr = nullptr;

            // Convert to SearchGraph weight pointers.
            const auto& [node_idx_0, weight_idx_0] = flooder_graph_weight_ptr_map.at(flooder_edge0_ptr);
            search_edge0_ptr = &(nodes[node_idx_0].neighbor_weights[weight_idx_0]);

            if (flooder_edge1_ptr != nullptr) {
                const auto& [node_idx_1, weight_idx_1] = flooder_graph_weight_ptr_map.at(flooder_edge1_ptr);
                search_edge1_ptr = &(nodes[node_idx_1].neighbor_weights[weight_idx_1]);
            }
            if (!previous_weights.count(search_edge0_ptr)) {
                previous_weights[search_edge0_ptr] = *search_edge0_ptr;
            }
            *search_edge0_ptr = std::min(*search_edge0_ptr, implied_weight);
            if (search_edge1_ptr != nullptr) {
                if (!previous_weights.count(search_edge1_ptr)) {
                    previous_weights[search_edge1_ptr] = *search_edge1_ptr;
                }
                *search_edge1_ptr = std::min(*search_edge1_ptr, implied_weight);
            }
        }
    }
}

void pm::SearchGraph::undo_reweights() {
    for (auto& [ptr, weight] : previous_weights) {
        *ptr = weight;
    }
    previous_weights.clear();
}
