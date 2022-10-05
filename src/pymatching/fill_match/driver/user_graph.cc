#include "pymatching/fill_match/driver/user_graph.h"

pm::UserNode::UserNode() {
}

void pm::UserGraph::add_edge(
    size_t node1, size_t node2, const std::vector<size_t>& observables, double weight, double error_probability) {
    if (node1 == node2)
        throw std::invalid_argument("Must have node1 != node2. Self loops not permitted.");
    auto max_id = std::max(node1, node2);
    if (max_id + 1 > nodes.size()) {
        nodes.reserve(2 * (max_id + 1));  // Ensure we don't allocate too often
        nodes.resize(max_id + 1);
    }
    nodes[node1].neighbors.push_back({node2, observables, weight, error_probability});
    nodes[node2].neighbors.push_back({node1, observables, weight, error_probability});

    for (auto& obs : observables) {
        if (obs + 1 > _num_observables)
            _num_observables = obs + 1;
    }
    _mwpm_needs_updating = true;

    if (error_probability < 0 || error_probability > 1)
        _all_edges_have_error_probabilities = false;
}

void pm::UserGraph::add_boundary_edge(
    size_t node, const std::vector<size_t>& observables, double weight, double error_probability) {
    if (node + 1 > nodes.size()) {
        nodes.reserve(2 * (node + 1));  // Ensure we don't allocate too often
        nodes.resize(node + 1);
    }
    nodes[node].neighbors.push_back({SIZE_MAX, observables, weight, error_probability});

    for (auto& obs : observables) {
        if (obs + 1 > _num_observables)
            _num_observables = obs + 1;
    }
    _mwpm_needs_updating = true;

    if (error_probability < 0 || error_probability > 1)
        _all_edges_have_error_probabilities = false;
}

pm::UserGraph::UserGraph()
    : _num_observables(0), _mwpm_needs_updating(true), _all_edges_have_error_probabilities(true) {
}

pm::UserGraph::UserGraph(size_t num_nodes)
    : _num_observables(0), _mwpm_needs_updating(true), _all_edges_have_error_probabilities(true) {
    nodes.resize(num_nodes);
}

pm::UserGraph::UserGraph(size_t num_nodes, size_t num_observables)
    : _num_observables(num_observables), _mwpm_needs_updating(true), _all_edges_have_error_probabilities(true) {
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

void pm::UserGraph::add_noise(uint8_t* error_arr, uint8_t* syndrome_arr) const {
    if (!_all_edges_have_error_probabilities)
        return;

    for (size_t i = 0; i < nodes.size(); i++) {
        for (size_t j = 0; j < nodes[i].neighbors.size(); j++) {
            size_t k = nodes[i].neighbors[j].node;
            if (i < k) {
                auto p = nodes[i].neighbors[j].error_probability;
                if (rand_float(0.0, 1.0) < p) {
                    // Flip the observables
                    for (auto& obs : nodes[i].neighbors[j].observable_indices) {
                        *(error_arr + obs) ^= 1;
                    }
                    // Flip the syndrome bits
                    *(syndrome_arr + i) ^= 1;
                    if (k != SIZE_MAX)
                        *(syndrome_arr + k) ^= 1;
                }
            }
        }
    }

    for (auto& b : boundary_nodes)
        *(syndrome_arr + b) = 0;
}
