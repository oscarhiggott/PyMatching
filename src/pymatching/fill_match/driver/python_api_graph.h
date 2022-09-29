#ifndef PYMATCHING2_PYTHON_API_GRAPH_H
#define PYMATCHING2_PYTHON_API_GRAPH_H

#include <set>
#include <vector>
#include <cmath>

#include "pymatching/fill_match/ints.h"

namespace pm {






struct UserNodeNeighbor {
    size_t node;                             /// The index in UserGraph.nodes of the neighboring node
    std::vector<size_t> observable_indices;  /// The indices of the observables crossed along this edge
    double neighbor_weight;                  /// The weight of the edge to this neighboring node
    double error_probability;                /// The error probability associated with this node
};

class UserNode {
   public:
    UserNode() {
    }
    std::vector<UserNodeNeighbor> neighbors;  /// The node's neighbors.
};

class UserGraph {
   public:
    std::vector<UserNode> nodes;
    std::set<size_t> boundary_nodes;

    UserGraph();
    void add_edge(size_t u, size_t v, const std::vector<size_t>& observables, double weight, double error_probability);
    void add_boundary_edge(size_t u, const std::vector<size_t>& observables, double weight, double error_probability);
};

}  // namespace pm

#endif  // PYMATCHING2_PYTHON_API_GRAPH_H
