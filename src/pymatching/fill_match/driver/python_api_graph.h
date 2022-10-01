#ifndef PYMATCHING2_PYTHON_API_GRAPH_H
#define PYMATCHING2_PYTHON_API_GRAPH_H

#include <cmath>
#include <set>
#include <vector>

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

class UserGraph {
   public:
    std::vector<UserNode> nodes;
    std::set<size_t> boundary_nodes;

    UserGraph();
    explicit UserGraph(size_t num_nodes);
    UserGraph(size_t num_nodes, size_t num_observables);
    void add_edge(size_t u, size_t v, const std::vector<size_t>& observables, double weight, double error_probability);
    void add_boundary_edge(size_t u, const std::vector<size_t>& observables, double weight, double error_probability);
    void set_boundary(const std::set<size_t>& boundary);
    std::set<size_t> get_boundary();
    size_t get_num_observables();
    size_t get_num_nodes();
    size_t get_num_detectors();
    bool is_boundary_node(size_t node_id);
    pm::IntermediateWeightedGraph to_intermediate_weighted_graph();
    void update_mwpm();
    Mwpm& get_mwpm();
   private:
    pm::Mwpm _mwpm;
    size_t _num_observables;
    bool _mwpm_needs_updating;
};

}  // namespace pm

#endif  // PYMATCHING2_PYTHON_API_GRAPH_H
