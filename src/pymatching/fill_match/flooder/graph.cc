#include "pymatching/fill_match/flooder/graph.h"

#include "pymatching/fill_match/flooder/graph_fill_region.h"
#include "pymatching/fill_match/flooder_matcher_interop/mwpm_event.h"

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

    pm::obs_int obs_mask = 0;
    if (num_observables <= sizeof(pm::obs_int) * 8){
        for (auto obs : observables)
            obs_mask ^= (pm::obs_int) 1 << obs;
    }

    if (weight < 0) {
        update_negative_weight_observables(observables);
        update_negative_weight_detection_events(u);
        update_negative_weight_detection_events(v);
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
            obs_mask ^= (pm::obs_int) 1 << obs;
    }

    if (weight < 0) {
        update_negative_weight_observables(observables);
        update_negative_weight_detection_events(u);
    }

    auto &n = nodes[u];
    if (!n.neighbors.empty() && n.neighbors[0] == nullptr) {
        throw std::invalid_argument("Max one boundary edge.");
    }
    n.neighbors.insert(n.neighbors.begin(), 1, nullptr);
    n.neighbor_weights.insert(n.neighbor_weights.begin(), 1, std::abs(weight));
    n.neighbor_observables.insert(n.neighbor_observables.begin(), 1, obs_mask);
}

MatchingGraph::MatchingGraph(size_t num_nodes, size_t num_observables)
    : num_nodes(num_nodes), num_observables(num_observables) {
    nodes.resize(num_nodes);
}

MatchingGraph::MatchingGraph(MatchingGraph &&graph) noexcept
    : nodes(std::move(graph.nodes)), num_nodes(graph.num_nodes), num_observables(graph.num_observables) {
}

MatchingGraph::MatchingGraph() : num_nodes(0), num_observables(0) {
}

void MatchingGraph::update_negative_weight_observables(const std::vector<size_t> &observables) {
    for (auto& obs : observables) {
        auto it = negative_weight_observables.find(obs);
        if (it == negative_weight_observables.end()) {
            negative_weight_observables.insert(obs);
        } else {
            negative_weight_observables.erase(it);
        }
    }
}

void MatchingGraph::update_negative_weight_detection_events(size_t node_id) {
    auto it = negative_weight_detection_events.find(node_id);
    if (it == negative_weight_detection_events.end()) {
        negative_weight_detection_events.insert(node_id);
    } else {
        negative_weight_detection_events.erase(it);
    }
}

}  // namespace pm
