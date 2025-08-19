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

#include <cassert>
#include <cmath>
#include <map>

#include "graph.h"
#include "pymatching/sparse_blossom/driver/implied_weights.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.h"

namespace pm {

void MatchingGraph::add_edge(
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

    // all_edges_to_implied_weights_unconverted[u][v] for a node u corresponds to the edge weights conditioned by (u, v)
    // where v is the v'th neighbour of u in nodes[u].neighbors.

    nodes[u].neighbors.push_back(&(nodes[v]));
    nodes[u].neighbor_weights.push_back(std::abs(weight));
    nodes[u].neighbor_observables.push_back(obs_mask);
    nodes[u].neighbor_implied_weights.push_back({});
    edges_to_implied_weights_unconverted[u].emplace_back(implied_weights_for_other_edges);

    nodes[v].neighbors.push_back(&(nodes[u]));
    nodes[v].neighbor_weights.push_back(std::abs(weight));
    nodes[v].neighbor_observables.push_back(obs_mask);
    nodes[v].neighbor_implied_weights.push_back({});
    edges_to_implied_weights_unconverted[v].emplace_back(implied_weights_for_other_edges);
}

void MatchingGraph::add_boundary_edge(
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
    n.neighbor_implied_weights.insert(n.neighbor_implied_weights.begin(), 1, {});
    edges_to_implied_weights_unconverted[u].insert(
        edges_to_implied_weights_unconverted[u].begin(), 1, implied_weights_for_other_edges);
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
      is_user_graph_boundary_node(std::move(graph.is_user_graph_boundary_node)),
      num_nodes(graph.num_nodes),
      num_observables(graph.num_observables),
      normalising_constant(graph.normalising_constant),
      previous_weights(graph.previous_weights),
      edges_to_implied_weights_unconverted(graph.edges_to_implied_weights_unconverted),
      loaded_from_dem_without_correlations(graph.loaded_from_dem_without_correlations) {
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

namespace {

ImpliedWeight convert_rule(
    std::vector<DetectorNode>& nodes, const ImpliedWeightUnconverted& rule, const double normalising_constant) {
    const size_t& i = rule.node1;
    const size_t& j = rule.node2;
    weight_int* weight_pointer_i =
        &nodes[i].neighbor_weights[nodes[i].index_of_neighbor(j == SIZE_MAX ? nullptr : &nodes[j])];
    weight_int* weight_pointer_j =
        j == SIZE_MAX ? nullptr : &nodes[j].neighbor_weights[nodes[j].index_of_neighbor(&nodes[i])];

    double rescaled_normalising_constant = normalising_constant / 2;
    pm::signed_weight_int w = (pm::signed_weight_int)round(rule.implied_weight * rescaled_normalising_constant);
    // Extremely important!
    // If all edge weights are even integers, then all collision events occur at integer times.
    w *= 2;

    return ImpliedWeight{weight_pointer_i, weight_pointer_j, static_cast<pm::weight_int>(std::abs(w))};
}

}  // namespace

void MatchingGraph::convert_implied_weights(double normalising_constant) {
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

// Reweight assuming an error has occurred on a single edge u, v. When v == -1, assumes an edge from
// u to the boundary.
void MatchingGraph::reweight_for_edge(const int64_t& u, const int64_t& v) {
    size_t z = nodes[u].index_of_neighbor(v == -1 ? nullptr : &nodes[v]);
    reweight(nodes[u].neighbor_implied_weights[z]);
}

void MatchingGraph::reweight_for_edges(const std::vector<int64_t>& edges) {
    if (loaded_from_dem_without_correlations) {
        throw std::invalid_argument(
            "Attempting to decode with `enable_correlations=True`, however the decoder has "
            "not been configured to decode using correlations. Ensure that you also set "
            "`enable_correlations=True` when loading the stim circuit or detector error model. "
            "For example: `matching = pymatching.Matching.from_detector_error_model(dem, enable_correlations=True)`");
    }
    for (size_t i = 0; i < edges.size() >> 1; ++i) {
        int64_t u = edges[2 * i];
        int64_t v = edges[2 * i + 1];
        reweight_for_edge(u, v);
    }
}

void MatchingGraph::undo_reweights() {
    // We iterate backward over the previous weights, since some edges
    // may have been reweighted more than once. Alternatively,
    // we could iterate forward and only undo a reweight if the
    // previous weight is larger.
    for (auto it = previous_weights.rbegin(); it != previous_weights.rend(); ++it) {
        pm::PreviousWeight& prev = *it;
        *prev.ptr = prev.val;
    }
    previous_weights.clear();
}

}  // namespace pm
