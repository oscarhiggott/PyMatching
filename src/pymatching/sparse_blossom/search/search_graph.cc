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

pm::SearchGraph::SearchGraph() : num_nodes(0) {
}

pm::SearchGraph::SearchGraph(size_t num_nodes) : num_nodes(num_nodes) {
    nodes.resize(num_nodes);
}

pm::SearchGraph::SearchGraph(pm::SearchGraph &&graph) noexcept
    : nodes(std::move(graph.nodes)), num_nodes(graph.num_nodes) {
}

void pm::SearchGraph::add_edge(size_t u, size_t v, signed_weight_int weight, const std::vector<size_t> &observables) {
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

    nodes[u].neighbors.push_back(&(nodes[v]));
    nodes[u].neighbor_weights.push_back(std::abs(weight));
    nodes[u].neighbor_observable_indices.push_back(observables);

    nodes[v].neighbors.push_back(&(nodes[u]));
    nodes[v].neighbor_weights.push_back(std::abs(weight));
    nodes[v].neighbor_observable_indices.push_back(observables);
}

void pm::SearchGraph::add_boundary_edge(size_t u, signed_weight_int weight, const std::vector<size_t> &observables) {
    if (u >= nodes.size()) {
        throw std::invalid_argument(
            "Node " + std::to_string(u) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }

    nodes[u].neighbors.insert(nodes[u].neighbors.begin(), 1, nullptr);
    nodes[u].neighbor_weights.insert(nodes[u].neighbor_weights.begin(), 1, std::abs(weight));
    nodes[u].neighbor_observable_indices.insert(nodes[u].neighbor_observable_indices.begin(), 1, observables);
}
