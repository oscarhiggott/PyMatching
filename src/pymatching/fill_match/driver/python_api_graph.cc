#include "pymatching/fill_match/driver/python_api_graph.h"

pm::UserNode::UserNode() {
}

void pm::UserGraph::add_edge(
    size_t u, size_t v, const std::vector<size_t>& observables, double weight, double error_probability) {
    if (u == v)
        throw std::invalid_argument("Must have u != v. Self loops not permitted.");
    auto max_id = std::max(u, v);
    if (max_id > nodes.size()) {
        nodes.reserve(2 * (max_id + 1));  // Ensure we don't allocate too often
        nodes.resize(max_id + 1);
    }
    nodes[u].neighbors.push_back({v, observables, weight, error_probability});
    nodes[v].neighbors.push_back({u, observables, weight, error_probability});

    for (auto& obs : observables) {
        if (obs + 1 > _num_observables)
            _num_observables = obs + 1;
    }
}

void pm::UserGraph::add_boundary_edge(
    size_t u, const std::vector<size_t>& observables, double weight, double error_probability) {
    if (u > nodes.size()) {
        nodes.reserve(2 * (u + 1));  // Ensure we don't allocate too often
        nodes.resize(u + 1);
    }
    nodes[u].neighbors.push_back({SIZE_MAX, observables, weight, error_probability});

    for (auto& obs : observables) {
        if (obs + 1 > _num_observables)
            _num_observables = obs + 1;
    }
}

pm::UserGraph::UserGraph() : _num_observables(0) {
}

pm::UserGraph::UserGraph(size_t num_nodes) : _num_observables(0) {
    nodes.resize(num_nodes);
}

pm::UserGraph::UserGraph(size_t num_nodes, size_t num_observables) : _num_observables(num_observables) {
    nodes.resize(num_nodes);
}

void pm::UserGraph::set_boundary(std::set<size_t>& boundary) {
    boundary_nodes = boundary;
}

std::set<size_t> pm::UserGraph::get_boundary() {
    return boundary_nodes;
}

pm::IntermediateWeightedGraph pm::UserGraph::to_intermediate_weighted_graph() {
    pm::IntermediateWeightedGraph graph(nodes.size(), get_num_observables());
    for (size_t i = 0; i < nodes.size(); i++) {
        for (size_t j = 0; j < nodes[i].neighbors.size(); j++) {
            auto n = nodes[i].neighbors[j];
            size_t v = n.node;
            if (i < v) {
                bool i_is_boundary = is_boundary_node(i);
                bool v_is_boundary = is_boundary_node(v);
                if (i_is_boundary && !v_is_boundary) {
                    graph.add_or_merge_boundary_edge(v, n.weight, n.observable_indices);
                } else if (v_is_boundary && !i_is_boundary) {
                    graph.add_or_merge_boundary_edge(i, n.weight, n.observable_indices);
                } else if (!v_is_boundary) {
                    // Neither i nor v are boundary nodes
                    graph.add_or_merge_edge(i, v, n.weight, n.observable_indices);
                }
            }
        }
    }
    return graph;
}

size_t pm::UserGraph::get_num_observables() {
    return _num_observables;
}

bool pm::UserGraph::is_boundary_node(size_t node_id) {
    return (node_id == SIZE_MAX) || (boundary_nodes.find(node_id) != boundary_nodes.end());
}
pm::Mwpm pm::UserGraph::to_mwpm() {
    auto graph = to_intermediate_weighted_graph();
    return graph.to_mwpm(1 << 14);
}
