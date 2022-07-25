#include<algorithm>
#include "graph.h"

namespace pm {

    template<template<class> class neighbor_vec>
    void Graph<neighbor_vec>::add_edge(size_t u, size_t v, weight_int weight, obs_int observables) {
        size_t larger_node = std::max(u, v);
        if (larger_node > nodes.size() - 1) {
            nodes.resize(larger_node);
        }
        // Allow parallel edges?
        nodes[u].neighbors.push_back(&nodes[v]);
        nodes[u].neighbor_weights.push_back(&nodes[v]);
        nodes[u].neighbor_observables.push_back(&nodes[v]);
        nodes[u].neighbor_schedules.push_back(nullptr);

        nodes[v].neighbors.push_back(&nodes[u]);
        nodes[v].neighbor_weights.push_back(&nodes[u]);
        nodes[v].neighbor_observables.push_back(&nodes[u]);
        nodes[v].neighbor_schedules.push_back(nullptr);
    }

}
