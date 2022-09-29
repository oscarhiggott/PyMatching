#ifndef PYMATCHING2_PYTHON_API_GRAPH_H
#define PYMATCHING2_PYTHON_API_GRAPH_H

#include <set>
#include <vector>
#include <cmath>

#include "pymatching/fill_match/ints.h"

namespace pm {





/// Computes the weight of an edge resulting from merging edges with weight `a' and weight `b', assuming each edge
/// weight is a log-likelihood ratio log((1-p)/p) associated with the probability p of an error occurring on the
/// edge, and that the error mechanisms associated with the two edges being merged are independent.
///
/// A mathematically equivalent implementation of this method would be:
///
///  double merge_weights(double a, double b){
///     double p_a = 1/(1 + std::exp(a));
///     double p_b = 1/(1 + std::exp(b));
///     double p_both = p_a * (1 - p_b) + p_b * (1 - p_a);
///     return std::log((1-p_both)/p_both);
///  }
///
/// however this would suffer from numerical overflow issues for abs(a) >> 30 or abs(b) >> 30 which is avoided by the
/// the implementation used here instead. See Equation (6) of https://ieeexplore.ieee.org/document/1495850 for more
/// details.
double merge_weights(double a, double b);

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
