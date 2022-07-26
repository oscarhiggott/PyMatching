#include<algorithm>
#include "graph.h"

namespace pm {

    void Graph::add_edge(size_t u, size_t v, weight_int weight, obs_int observables) {
        size_t larger_node = std::max(u, v);
        if (larger_node + 1 > nodes.size())
            nodes.resize(larger_node + 1);

        // Allow parallel edges?
        nodes[u].neighbors.push_back(&(nodes[v]));
        nodes[u].neighbor_weights.push_back(weight);
        nodes[u].neighbor_observables.push_back(observables);
        nodes[u].neighbor_schedules.push_back(nullptr);

        nodes[v].neighbors.push_back(&(nodes[u]));
        nodes[v].neighbor_weights.push_back(weight);
        nodes[v].neighbor_observables.push_back(observables);
        nodes[v].neighbor_schedules.push_back(nullptr);
    }

    void Graph::add_boundary_edge(size_t u, weight_int weight, obs_int observables) {
        if (u > nodes.size() - 1)
            nodes.resize(u + 1);

        nodes[u].neighbors.push_back(nullptr);
        nodes[u].neighbor_weights.push_back(weight);
        nodes[u].neighbor_observables.push_back(observables);
        nodes[u].neighbor_schedules.push_back(nullptr);
    }

    Graph::Graph(size_t num_nodes) : num_nodes(num_nodes) {
        nodes.resize(num_nodes);
    }

}
