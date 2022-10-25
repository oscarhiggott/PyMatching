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

#include "pymatching/sparse_blossom/flooder/graph.h"

#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.h"

namespace pm {

void MatchingGraph::add_edge(size_t u, size_t v, signed_weight_int weight, const std::vector<size_t>& observables) {
    size_t larger_node = std::max(u, v);
    if (larger_node + 1 > nodes.size()) {
        throw std::invalid_argument(
            "Node " + std::to_string(larger_node) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }

    if (weight < 0) {
        update_negative_weight_observables(observables);
        update_negative_weight_detection_events(u);
        update_negative_weight_detection_events(v);
        negative_weight_sum += weight;
    }

    // We don't actually add self-loops into the graph (since a self loop with positive weight should never be added
    // to the solution, since it causes no detection events). However, if it has a negative edge weight, then it is
    // always included, and is handled by the negative weight case above.
    if (u == v)
        return;

    pm::obs_int obs_mask = 0;
    if (num_observables <= sizeof(pm::obs_int) * 8) {
        for (auto obs : observables)
            obs_mask ^= (pm::obs_int)1 << obs;
    }

    nodes[u].neighbors.push_back(&(nodes[v]));
    nodes[u].neighbor_weights.push_back(std::abs(weight));
    nodes[u].neighbor_observables.push_back(obs_mask);

    nodes[v].neighbors.push_back(&(nodes[u]));
    nodes[v].neighbor_weights.push_back(std::abs(weight));
    nodes[v].neighbor_observables.push_back(obs_mask);
}

void MatchingGraph::add_boundary_edge(size_t u, signed_weight_int weight, const std::vector<size_t>& observables) {
    if (u >= nodes.size()) {
        throw std::invalid_argument(
            "Node " + std::to_string(u) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }

    pm::obs_int obs_mask = 0;
    if (num_observables <= sizeof(pm::obs_int) * 8) {
        for (auto obs : observables)
            obs_mask ^= (pm::obs_int)1 << obs;
    }

    if (weight < 0) {
        update_negative_weight_observables(observables);
        update_negative_weight_detection_events(u);
        negative_weight_sum += weight;
    }

    auto& n = nodes[u];
    if (!n.neighbors.empty() && n.neighbors[0] == nullptr) {
        throw std::invalid_argument("Max one boundary edge.");
    }
    n.neighbors.insert(n.neighbors.begin(), 1, nullptr);
    n.neighbor_weights.insert(n.neighbor_weights.begin(), 1, std::abs(weight));
    n.neighbor_observables.insert(n.neighbor_observables.begin(), 1, obs_mask);
}

MatchingGraph::MatchingGraph(size_t num_nodes, size_t num_observables)
    : negative_weight_sum(0), num_nodes(num_nodes), num_observables(num_observables), normalising_constant(0) {
    nodes.resize(num_nodes);
}

MatchingGraph::MatchingGraph(size_t num_nodes, size_t num_observables, double normalising_constant)
    : negative_weight_sum(0),
      num_nodes(num_nodes),
      num_observables(num_observables),
      normalising_constant(normalising_constant) {
    nodes.resize(num_nodes);
}

MatchingGraph::MatchingGraph(MatchingGraph&& graph) noexcept
    : nodes(std::move(graph.nodes)),
      negative_weight_detection_events_set(std::move(graph.negative_weight_detection_events_set)),
      negative_weight_observables_set(std::move(graph.negative_weight_observables_set)),
      negative_weight_sum(graph.negative_weight_sum),
      num_nodes(graph.num_nodes),
      num_observables(graph.num_observables),
      normalising_constant(graph.normalising_constant) {
}

MatchingGraph::MatchingGraph() : negative_weight_sum(0), num_nodes(0), num_observables(0), normalising_constant(0) {
}

void MatchingGraph::update_negative_weight_observables(const std::vector<size_t>& observables) {
    for (auto& obs : observables) {
        auto it = negative_weight_observables_set.find(obs);
        if (it == negative_weight_observables_set.end()) {
            negative_weight_observables_set.insert(obs);
        } else {
            negative_weight_observables_set.erase(it);
        }
    }
}

void MatchingGraph::update_negative_weight_detection_events(size_t node_id) {
    auto it = negative_weight_detection_events_set.find(node_id);
    if (it == negative_weight_detection_events_set.end()) {
        negative_weight_detection_events_set.insert(node_id);
    } else {
        negative_weight_detection_events_set.erase(it);
    }
}

}  // namespace pm
