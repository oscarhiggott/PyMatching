#ifndef PYMATCHING2_USER_GRAPH_H
#define PYMATCHING2_USER_GRAPH_H

#include <cmath>
#include <set>
#include <vector>
#include "pymatching/rand/rand_gen.h"

#include "pymatching/fill_match/driver/io.h"
#include "pymatching/fill_match/ints.h"


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
    std::vector<UserNodeNeighbor> neighbors;  /// The node's neighbors.
};

typedef std::tuple<size_t,size_t,std::vector<size_t>,double, double> edge_data;

class UserGraph {
   public:
    std::vector<UserNode> nodes;
    std::set<size_t> boundary_nodes;

    UserGraph();
    explicit UserGraph(size_t num_nodes);
    UserGraph(size_t num_nodes, size_t num_observables);
    void add_edge(size_t node1, size_t node2, const std::vector<size_t>& observables, double weight, double error_probability);
    void add_boundary_edge(size_t node, const std::vector<size_t>& observables, double weight, double error_probability);
    void set_boundary(const std::set<size_t>& boundary);
    std::set<size_t> get_boundary();
    size_t get_num_observables();
    size_t get_num_nodes();
    size_t get_num_detectors();
    size_t get_num_edges();
    bool is_boundary_node(size_t node_id);
    pm::IntermediateWeightedGraph to_intermediate_weighted_graph();
    void update_mwpm();
    Mwpm& get_mwpm();
    void add_noise(uint8_t* error_arr, uint8_t* syndrome_arr) const;
    bool all_edges_have_error_probabilities();
    std::vector<edge_data> get_edges();
   private:
    pm::Mwpm _mwpm;
    size_t _num_observables;
    size_t _num_edges;
    bool _mwpm_needs_updating;
    bool _all_edges_have_error_probabilities;
};

}  // namespace pm

#endif  // PYMATCHING2_USER_GRAPH_H
