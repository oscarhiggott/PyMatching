#include "pymatching/fill_match/flooder/graph.h"

#include "pymatching/fill_match/flooder/graph_fill_region.h"
#include "pymatching/fill_match/flooder_matcher_interop/mwpm_event.h"

namespace pm {

void MatchingGraph::add_edge(size_t u, size_t v, weight_int weight, obs_int observables) {
    size_t larger_node = std::max(u, v);
    if (larger_node + 1 > nodes.size()) {
        throw std::invalid_argument(
            "Node " + std::to_string(larger_node) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }

    // Allow parallel edges?
    nodes[u].neighbors.push_back(&(nodes[v]));
    nodes[u].neighbor_weights.push_back(weight);
    nodes[u].neighbor_observables.push_back(observables);

    nodes[v].neighbors.push_back(&(nodes[u]));
    nodes[v].neighbor_weights.push_back(weight);
    nodes[v].neighbor_observables.push_back(observables);
}

void MatchingGraph::add_boundary_edge(size_t u, weight_int weight, obs_int observables) {
    if (u >= nodes.size()) {
        throw std::invalid_argument(
            "Node " + std::to_string(u) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }
    auto &n = nodes[u];
    if (!n.neighbors.empty() && n.neighbors[0] == nullptr) {
        throw std::invalid_argument("Max one boundary edge.");
    }
    n.neighbors.insert(n.neighbors.begin(), 1, nullptr);
    n.neighbor_weights.insert(n.neighbor_weights.begin(), 1, weight);
    n.neighbor_observables.insert(n.neighbor_observables.begin(), 1, observables);
}

MatchingGraph::MatchingGraph(size_t num_nodes) : num_nodes(num_nodes) {
    nodes.resize(num_nodes);
}

MatchingGraph::MatchingGraph(MatchingGraph &&graph) noexcept
    : nodes(std::move(graph.nodes)), num_nodes(graph.num_nodes) {
}

MatchingGraph::MatchingGraph() : num_nodes(0) {
}

}  // namespace pm
