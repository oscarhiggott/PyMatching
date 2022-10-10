#include "pymatching/fill_match/driver/user_graph.h"

pm::UserNode::UserNode() : is_boundary(false) {
}

size_t pm::UserNode::index_of_neighbor(size_t node) const {
    auto it = std::find_if(neighbors.begin(), neighbors.end(), [&](const UserNodeNeighbor& neighbor) {
        return neighbor.node == node;
    });
    if (it == neighbors.end())
        return SIZE_MAX;

    return it - neighbors.begin();
}

bool is_valid_probability(double p) {
    return (p >= 0 && p <= 1);
}

void pm::UserGraph::merge_edge_or_boundary_edge(
    size_t node,
    size_t neighbor_index,
    const std::vector<size_t>& parallel_observables,
    double parallel_weight,
    double parallel_error_probability,
    pm::MERGE_STRATEGY merge_strategy) {
    auto& neighbor = nodes[node].neighbors[neighbor_index];
    if (merge_strategy == DISALLOW) {
        throw std::invalid_argument("Parallel edges not permitted with provided merge strategy");
    } else if (
        merge_strategy == KEEP_ORIGINAL || (merge_strategy == SMALLEST_WEIGHT && parallel_weight >= neighbor.weight)) {
        return;
    } else {
        double new_weight, new_error_probability;
        bool use_new_observables;
        if (merge_strategy == REPLACE || merge_strategy == SMALLEST_WEIGHT) {
            new_weight = parallel_weight;
            new_error_probability = parallel_error_probability;
            use_new_observables = true;
        } else if (merge_strategy == INDEPENDENT) {
            new_weight = pm::merge_weights(parallel_weight, neighbor.weight);
            new_error_probability = -1;
            if (is_valid_probability(neighbor.error_probability) && is_valid_probability(parallel_error_probability))
                new_error_probability = parallel_error_probability * (1 - neighbor.error_probability) +
                                        neighbor.error_probability * (1 - parallel_error_probability);
            // We do not need to update the observables. If they do not match up, then the code has distance 2.
            use_new_observables = false;
        } else {
            throw std::invalid_argument("Merge strategy not recognised.");
        }
        // Update the existing edge weight and probability in the adjacency list of `node`
        neighbor.weight = new_weight;
        neighbor.error_probability = new_error_probability;
        if (use_new_observables)
            neighbor.observable_indices = parallel_observables;

        if (neighbor.node != SIZE_MAX) {
            // Update the existing edge weight and probability in the adjacency list of node2 if not the boundary
            size_t idx2 = nodes[neighbor.node].index_of_neighbor(node);
            nodes[neighbor.node].neighbors[idx2].weight = new_weight;
            nodes[neighbor.node].neighbors[idx2].error_probability = new_error_probability;
            if (use_new_observables)
                nodes[neighbor.node].neighbors[idx2].observable_indices = parallel_observables;
        }

        _mwpm_needs_updating = true;
        if (new_error_probability < 0 || new_error_probability > 1)
            _all_edges_have_error_probabilities = false;
    }
}

void pm::UserGraph::add_or_merge_edge(
    size_t node1,
    size_t node2,
    const std::vector<size_t>& observables,
    double weight,
    double error_probability,
    MERGE_STRATEGY merge_strategy) {
    if (node1 == node2)
        throw std::invalid_argument("Must have node1 != node2. Self loops not permitted.");
    auto max_id = std::max(node1, node2);
    if (max_id + 1 > nodes.size())
        nodes.resize(max_id + 1);

    size_t idx = nodes[node1].index_of_neighbor(node2);

    if (idx == SIZE_MAX) {
        nodes[node1].neighbors.push_back({node2, observables, weight, error_probability});
        nodes[node2].neighbors.push_back({node1, observables, weight, error_probability});

        for (auto& obs : observables) {
            if (obs + 1 > _num_observables)
                _num_observables = obs + 1;
        }
        _num_edges++;
        _mwpm_needs_updating = true;
        if (error_probability < 0 || error_probability > 1)
            _all_edges_have_error_probabilities = false;
    } else {
        merge_edge_or_boundary_edge(node1, idx, observables, weight, error_probability, merge_strategy);
    }
}

void pm::UserGraph::add_or_merge_boundary_edge(
    size_t node,
    const std::vector<size_t>& observables,
    double weight,
    double error_probability,
    MERGE_STRATEGY merge_strategy) {
    if (node + 1 > nodes.size())
        nodes.resize(node + 1);

    size_t idx = nodes[node].index_of_neighbor(SIZE_MAX);

    if (idx == SIZE_MAX) {
        nodes[node].neighbors.push_back({SIZE_MAX, observables, weight, error_probability});

        for (auto& obs : observables) {
            if (obs + 1 > _num_observables)
                _num_observables = obs + 1;
        }
        _num_edges++;
        _has_pseudo_boundary_node = true;
        _mwpm_needs_updating = true;
        if (error_probability < 0 || error_probability > 1)
            _all_edges_have_error_probabilities = false;
    } else {
        merge_edge_or_boundary_edge(node, idx, observables, weight, error_probability, merge_strategy);
    }
}

pm::UserGraph::UserGraph()
    : _num_observables(0),
      _num_edges(0),
      _mwpm_needs_updating(true),
      _all_edges_have_error_probabilities(true),
      _has_pseudo_boundary_node(false) {
}

pm::UserGraph::UserGraph(size_t num_nodes)
    : _num_observables(0),
      _num_edges(0),
      _mwpm_needs_updating(true),
      _all_edges_have_error_probabilities(true),
      _has_pseudo_boundary_node(false) {
    nodes.resize(num_nodes);
}

pm::UserGraph::UserGraph(size_t num_nodes, size_t num_observables)
    : _num_observables(num_observables),
      _num_edges(0),
      _mwpm_needs_updating(true),
      _all_edges_have_error_probabilities(true),
      _has_pseudo_boundary_node(false) {
    nodes.resize(num_nodes);
}

void pm::UserGraph::set_boundary(const std::set<size_t>& boundary) {
    for (auto& n : boundary_nodes)
        nodes[n].is_boundary = false;
    boundary_nodes = boundary;
    for (auto& n : boundary_nodes) {
        if (n >= nodes.size())
            nodes.resize(n + 1);
        nodes[n].is_boundary = true;
    }
    _mwpm_needs_updating = true;
}

std::set<size_t> pm::UserGraph::get_boundary() {
    return boundary_nodes;
}

size_t pm::UserGraph::get_num_observables() {
    return _num_observables;
}

size_t pm::UserGraph::get_num_nodes() {
    return nodes.size();
}

size_t pm::UserGraph::get_num_detectors() {
    return get_num_nodes() - boundary_nodes.size();
}

bool pm::UserGraph::is_boundary_node(size_t node_id) {
    return (node_id == SIZE_MAX) || nodes[node_id].is_boundary;
}

void pm::UserGraph::update_mwpm() {
    _mwpm = to_mwpm(pm::NUM_DISTINCT_WEIGHTS_FROM_USER_GRAPH, false);
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

size_t pm::UserGraph::get_num_edges() {
    return _num_edges;
}

bool pm::UserGraph::all_edges_have_error_probabilities() {
    return _all_edges_have_error_probabilities;
}
std::vector<pm::edge_data> pm::UserGraph::get_edges() {
    std::vector<pm::edge_data> edges;
    for (size_t i = 0; i < nodes.size(); i++) {
        for (size_t j = 0; j < nodes[i].neighbors.size(); j++) {
            auto n = nodes[i].neighbors[j];
            size_t v = n.node;
            if (i < v) {
                double p;
                if (n.error_probability < 0 || n.error_probability > 1) {
                    p = -1.0;
                } else {
                    p = n.error_probability;
                }
                if (v == SIZE_MAX)
                    v = nodes.size();  // Virtual boundary node has an index one larger than the largest current index
                edges.push_back({i, v, n.observable_indices, n.weight, p});
            }
        }
    }
    return edges;
}

double pm::UserGraph::max_abs_weight() {
    double max_abs_weight = 0;
    for (auto& node : nodes) {
        for (auto& neighbor : node.neighbors) {
            if (std::abs(neighbor.weight) > max_abs_weight)
                max_abs_weight = std::abs(neighbor.weight);
        }
    }
    return max_abs_weight;
}

pm::MatchingGraph pm::UserGraph::to_matching_graph(pm::weight_int num_distinct_weights) {
    pm::MatchingGraph matching_graph(nodes.size(), _num_observables);
    double normalising_constant = iter_discretized_edges(
        num_distinct_weights,
        [&](size_t u, size_t v, pm::signed_weight_int weight, const std::vector<size_t>& observables) {
            matching_graph.add_edge(u, v, weight, observables);
        },
        [&](size_t u, pm::signed_weight_int weight, const std::vector<size_t>& observables) {
            // Only add the boundary edge if it already isn't present. Ideally parallel edges should already have been
            // merged, however we are implicitly merging all boundary nodes in this step, which could give rise to new
            // parallel edges.
            if (matching_graph.nodes[u].neighbors.empty() || matching_graph.nodes[u].neighbors[0])
                matching_graph.add_boundary_edge(u, weight, observables);
        });
    matching_graph.normalising_constant = normalising_constant;
    return matching_graph;
}

pm::SearchGraph pm::UserGraph::to_search_graph(pm::weight_int num_distinct_weights) {
    /// Identical to to_matching_graph but for constructing a pm::SearchGraph
    pm::SearchGraph search_graph(nodes.size());
    iter_discretized_edges(
        num_distinct_weights,
        [&](size_t u, size_t v, pm::signed_weight_int weight, const std::vector<size_t>& observables) {
            search_graph.add_edge(u, v, weight, observables);
        },
        [&](size_t u, pm::signed_weight_int weight, const std::vector<size_t>& observables) {
            // Only add the boundary edge if it already isn't present. Ideally parallel edges should already have been
            // merged, however we are implicitly merging all boundary nodes in this step, which could give rise to new
            // parallel edges.
            if (search_graph.nodes[u].neighbors.empty() || search_graph.nodes[u].neighbors[0])
                search_graph.add_boundary_edge(u, weight, observables);
        });
    return search_graph;
}

pm::Mwpm pm::UserGraph::to_mwpm(pm::weight_int num_distinct_weights, bool ensure_search_graph_included) {
    if (_num_observables > sizeof(pm::obs_int) * 8 || ensure_search_graph_included) {
        auto mwpm = pm::Mwpm(
            pm::GraphFlooder(to_matching_graph(num_distinct_weights)),
            pm::SearchFlooder(to_search_graph(num_distinct_weights)));
        mwpm.flooder.sync_negative_weight_observables_and_detection_events();
        return mwpm;
    } else {
        auto mwpm = pm::Mwpm(pm::GraphFlooder(to_matching_graph(num_distinct_weights)));
        mwpm.flooder.sync_negative_weight_observables_and_detection_events();
        return mwpm;
    }
}

pm::Mwpm& pm::UserGraph::get_mwpm_with_search_graph() {
    if (!_mwpm_needs_updating && _mwpm.flooder.graph.nodes.size() == _mwpm.search_flooder.graph.nodes.size()) {
        return _mwpm;
    } else {
        _mwpm = to_mwpm(pm::NUM_DISTINCT_WEIGHTS_FROM_USER_GRAPH, true);
        _mwpm_needs_updating = false;
        return _mwpm;
    }
}

void pm::UserGraph::handle_dem_instruction(
    double p, const std::vector<size_t>& detectors, const std::vector<size_t>& observables) {
    if (detectors.size() == 2) {
        add_or_merge_edge(detectors[0], detectors[1], observables, std::log((1 - p) / p), p, INDEPENDENT);
    } else if (detectors.size() == 1) {
        add_or_merge_boundary_edge(detectors[0], observables, std::log((1 - p) / p), p, INDEPENDENT);
    }
}

void pm::UserGraph::get_nodes_on_shortest_path_from_source(size_t src, size_t dst, std::vector<size_t>& out_nodes) {
    auto& mwpm = get_mwpm_with_search_graph();
    bool src_is_boundary = is_boundary_node(src);
    bool dst_is_boundary = is_boundary_node(dst);
    if (src != SIZE_MAX && src >= nodes.size())
        throw std::invalid_argument("node " + std::to_string(src) + " is not in the graph");
    if (dst != SIZE_MAX && dst >= nodes.size())
        throw std::invalid_argument("node " + std::to_string(dst) + " is not in the graph");
    if (!src_is_boundary) {
        size_t dst_tmp = dst_is_boundary ? SIZE_MAX : dst;
        mwpm.search_flooder.iter_edges_on_shortest_path_from_source(src, dst_tmp, [&](const pm::SearchGraphEdge edge) {
            out_nodes.push_back(edge.detector_node - &mwpm.search_flooder.graph.nodes[0]);
        });
        if (!dst_is_boundary)
            out_nodes.push_back(dst);
    } else if (!dst_is_boundary) {
        std::vector<size_t> temp_out_nodes;
        get_nodes_on_shortest_path_from_source(dst, src, temp_out_nodes);
        for (size_t i = 0; i < temp_out_nodes.size(); i++) {
            out_nodes.push_back(temp_out_nodes[temp_out_nodes.size() - 1 - i]);
        }
    } else {
        throw std::invalid_argument("Both the source and destination vertices provided are boundary nodes");
    }
}

bool pm::UserGraph::has_edge(size_t node1, size_t node2) {
    if (node1 >= nodes.size())
        return false;
    return nodes[node1].index_of_neighbor(node2) != SIZE_MAX;
}

bool pm::UserGraph::has_boundary_edge(size_t node) {
    if (node >= nodes.size())
        return false;
    return nodes[node].index_of_neighbor(SIZE_MAX) != SIZE_MAX;
}

pm::UserGraph pm::detector_error_model_to_user_graph(const stim::DetectorErrorModel& detector_error_model) {
    pm::UserGraph user_graph(detector_error_model.count_detectors(), detector_error_model.count_observables());
    pm::iter_detector_error_model_edges(
        detector_error_model, [&](double p, const std::vector<size_t>& detectors, std::vector<size_t>& observables) {
            user_graph.handle_dem_instruction(p, detectors, observables);
        });
    return user_graph;
}