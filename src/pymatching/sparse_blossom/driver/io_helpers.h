#ifndef PYMATCHING2_IO_HELPERS_H_
#define PYMATCHING2_IO_HELPERS_H_

#include <cstddef>
#include <pymatching/sparse_blossom/ints.h>
#include <vector>

#include "pymatching/sparse_blossom/flooder/detector_node.h"
#include "stim.h"

namespace pm {

struct DetectorEdgeNodes {
    size_t d1;
    size_t d2;

    DetectorEdgeNodes();
    DetectorEdgeNodes(const size_t d);
    DetectorEdgeNodes(const size_t d1, const size_t d2);

    bool operator==(const DetectorEdgeNodes &other) const;
    bool operator!=(const DetectorEdgeNodes &other) const;
    bool operator<(const DetectorEdgeNodes &other) const;

    bool is_boundary() const;
};

struct UserEdge {
    DetectorEdgeNodes detectors;
    std::vector<size_t> observables;     /// The bservables crossed along this edge
    double weight;                       /// The weight of the edge to this neighboring node
    double error_probability;            /// The error probability associated with this node
    double correlated_proabilities_sum;  /// The error probability associated with this node
    std::vector<ImpliedWeightUnconverted> implied_weights_for_other_edges{};
    bool operator==(const UserEdge &other) const;
    bool operator!=(const UserEdge &other) const;
};

struct DecomposedDemError {
    /// The probability of this error occurring.
    double probability;
    /// Effects of the error. In normal surface code we never decompose into more
    /// than 4, but here it's 8 as a safety factor. It's generally not a great
    /// sign to see so many components, but it's more frustrating to have the
    /// decoding infrastructure refuse to even accept a wonky DEM.
    stim::FixedCapVector<UserEdge, 8> components;

    bool operator==(const DecomposedDemError &other) const;
    bool operator!=(const DecomposedDemError &other) const;
};

double bernoulli_xor_tmp(double p1, double p2);

void add_decomposed_error_to_conditional_groups(
    pm::DecomposedDemError &error,
    std::map<DetectorEdgeNodes, std::map<DetectorEdgeNodes, double>> &conditional_groups);

}  // namespace pm
#endif
