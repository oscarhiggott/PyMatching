#ifndef PYMATCHING2_USER_GRAPH_H
#define PYMATCHING2_USER_GRAPH_H

#include <cmath>
#include <set>
#include <vector>

#include "pymatching/fill_match/driver/io.h"
#include "pymatching/fill_match/ints.h"
#include "pymatching/rand/rand_gen.h"

namespace pm {

struct UserNodeNeighbor {
    size_t node;                             /// The index in UserGraph.nodes of the neighboring node
    std::vector<size_t> observable_indices;  /// The indices of the observables crossed along this edge
    double weight;                           /// The weight of the edge to this neighboring node
    double error_probability;                /// The error probability associated with this node
};

class UserNode {
   public:
    UserNode();
    size_t index_of_neighbor(size_t node) const;
    std::vector<UserNodeNeighbor> neighbors;  /// The node's neighbors.
    bool is_boundary;
};

typedef std::tuple<size_t, size_t, std::vector<size_t>, double, double> edge_data;

const pm::weight_int NUM_DISTINCT_WEIGHTS_FROM_USER_GRAPH = 1 << 14;

enum MERGE_STRATEGY : uint8_t { DISALLOW, INDEPENDENT, SMALLEST_WEIGHT, KEEP_ORIGINAL, REPLACE };

class UserGraph {
   public:
    std::vector<UserNode> nodes;
    std::set<size_t> boundary_nodes;

    UserGraph();
    explicit UserGraph(size_t num_nodes);
    UserGraph(size_t num_nodes, size_t num_observables);
    void merge_edge_or_boundary_edge(
        size_t node,
        size_t neighbor_index,
        const std::vector<size_t>& parallel_observables,
        double parallel_weight,
        double parallel_error_probability,
        MERGE_STRATEGY merge_strategy = DISALLOW);
    void add_or_merge_edge(
        size_t node1,
        size_t node2,
        const std::vector<size_t>& observables,
        double weight,
        double error_probability,
        MERGE_STRATEGY merge_strategy = DISALLOW);
    void add_or_merge_boundary_edge(
        size_t node,
        const std::vector<size_t>& observables,
        double weight,
        double error_probability,
        MERGE_STRATEGY merge_strategy = DISALLOW);
    bool has_edge(size_t node1, size_t node2);
    bool has_boundary_edge(size_t node);
    void set_boundary(const std::set<size_t>& boundary);
    std::set<size_t> get_boundary();
    size_t get_num_observables();
    size_t get_num_nodes();
    size_t get_num_detectors();
    size_t get_num_edges();
    bool is_boundary_node(size_t node_id);
    void add_noise(uint8_t* error_arr, uint8_t* syndrome_arr) const;
    bool all_edges_have_error_probabilities();
    std::vector<edge_data> get_edges();
    double max_abs_weight();
    template <typename EdgeCallable, typename BoundaryEdgeCallable>
    double iter_discretized_edges(
        pm::weight_int num_distinct_weights,
        const EdgeCallable& edge_func,
        const BoundaryEdgeCallable& boundary_edge_func);
    pm::MatchingGraph to_matching_graph(pm::weight_int num_distinct_weights);
    pm::SearchGraph to_search_graph(pm::weight_int num_distinct_weights);
    pm::Mwpm to_mwpm(pm::weight_int num_distinct_weights, bool ensure_search_graph_included);
    void update_mwpm();
    Mwpm& get_mwpm();
    Mwpm& get_mwpm_with_search_graph();
    void handle_dem_instruction(double p, const std::vector<size_t>& detectors, const std::vector<size_t>& observables);
    void get_nodes_on_shortest_path_from_source(size_t src, size_t dst, std::vector<size_t>& out_nodes);

   private:
    pm::Mwpm _mwpm;
    size_t _num_observables;
    size_t _num_edges;
    bool _mwpm_needs_updating;
    bool _all_edges_have_error_probabilities;
    bool _has_pseudo_boundary_node;
};

template <typename EdgeCallable, typename BoundaryEdgeCallable>
inline double UserGraph::iter_discretized_edges(
    pm::weight_int num_distinct_weights,
    const EdgeCallable& edge_func,
    const BoundaryEdgeCallable& boundary_edge_func) {
    double max_weight = max_abs_weight();
    pm::MatchingGraph matching_graph(nodes.size(), _num_observables);
    pm::weight_int max_half_edge_weight = num_distinct_weights - 1;
    double normalising_constant = (double)max_half_edge_weight / max_weight;
    for (auto& node : nodes) {
        for (auto& neighbor : node.neighbors) {
            auto i = &node - &nodes[0];
            pm::signed_weight_int w = (pm::signed_weight_int)(neighbor.weight * normalising_constant);

            // Extremely important!
            // If all edge weights are even integers, then all collision events occur at integer times.
            w *= 2;
            if (i < neighbor.node) {
                bool i_boundary = is_boundary_node(i);
                bool neighbor_boundary = is_boundary_node(neighbor.node);
                if (neighbor_boundary && !i_boundary) {
                    boundary_edge_func(i, w, neighbor.observable_indices);
                } else if (i_boundary && !neighbor_boundary) {
                    boundary_edge_func(neighbor.node, w, neighbor.observable_indices);
                } else if (!i_boundary) {
                    edge_func(i, neighbor.node, w, neighbor.observable_indices);
                }
            }
        }
    }
    return normalising_constant * 2;
}

UserGraph detector_error_model_to_user_graph(const stim::DetectorErrorModel& detector_error_model);

}  // namespace pm

#endif  // PYMATCHING2_USER_GRAPH_H
