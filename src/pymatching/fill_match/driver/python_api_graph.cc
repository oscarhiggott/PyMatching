#include "pymatching/fill_match/driver/python_api_graph.h"

pm::UserNode::UserNode() {
}

void pm::UserGraph::add_edge(
    size_t u, size_t v, const std::vector<size_t>& observables, double weight, double error_probability) {
    if (u == v)
        throw std::invalid_argument("Must have u != v. Self loops not permitted.");
    auto max_id = std::max(u, v);
    if (max_id + 1 > nodes.size()) {
        nodes.reserve(2 * (max_id + 1));  // Ensure we don't allocate too often
        nodes.resize(max_id + 1);
    }
    nodes[u].neighbors.push_back({v, observables, weight, error_probability});
    nodes[v].neighbors.push_back({u, observables, weight, error_probability});

    for (auto& obs : observables) {
        if (obs + 1 > _num_observables)
            _num_observables = obs + 1;
    }
    _mwpm_needs_updating = true;
}

void pm::UserGraph::add_boundary_edge(
    size_t u, const std::vector<size_t>& observables, double weight, double error_probability) {
    if (u + 1 > nodes.size()) {
        nodes.reserve(2 * (u + 1));  // Ensure we don't allocate too often
        nodes.resize(u + 1);
    }
    nodes[u].neighbors.push_back({SIZE_MAX, observables, weight, error_probability});

    for (auto& obs : observables) {
        if (obs + 1 > _num_observables)
            _num_observables = obs + 1;
    }
    _mwpm_needs_updating = true;
}

pm::UserGraph::UserGraph() : _num_observables(0), _mwpm_needs_updating(true) {
}

pm::UserGraph::UserGraph(size_t num_nodes) : _num_observables(0), _mwpm_needs_updating(true) {
    nodes.resize(num_nodes);
}

pm::UserGraph::UserGraph(size_t num_nodes, size_t num_observables)
    : _num_observables(num_observables), _mwpm_needs_updating(true) {
    nodes.resize(num_nodes);
}

void pm::UserGraph::set_boundary(const std::set<size_t>& boundary) {
    boundary_nodes = boundary;
    _mwpm_needs_updating = true;
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

size_t pm::UserGraph::get_num_nodes() {
    return nodes.size();
}

size_t pm::UserGraph::get_num_detectors() {
    return 0;
}

bool pm::UserGraph::is_boundary_node(size_t node_id) {
    return (node_id == SIZE_MAX) || (boundary_nodes.find(node_id) != boundary_nodes.end());
}

void pm::UserGraph::update_mwpm() {
    auto graph = to_intermediate_weighted_graph();
    _mwpm = std::move(graph.to_mwpm(1 << 14));
    _mwpm_needs_updating = false;
}
pm::Mwpm& pm::UserGraph::get_mwpm() {
    if (_mwpm_needs_updating)
        update_mwpm();
    return _mwpm;
}
