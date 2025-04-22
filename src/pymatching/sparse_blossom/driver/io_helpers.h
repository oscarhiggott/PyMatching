#ifndef PYMATCHING2_IO_HELPERS_H_
#define PYMATCHING2_IO_HELPERS_H_

#include <cstddef>
#include <pymatching/sparse_blossom/ints.h>
#include <vector>

#include "pymatching/sparse_blossom/flooder/detector_node.h"
#include "stim.h"

namespace pm {

struct DetectorEdgeId {
    size_t d1;
    size_t d2;

    DetectorEdgeId();
    DetectorEdgeId(const size_t d);
    DetectorEdgeId(const size_t d1, const size_t d2);

    bool operator==(const DetectorEdgeId &other) const;
    bool operator!=(const DetectorEdgeId &other) const;
    bool operator<(const DetectorEdgeId &other) const;

    bool is_boundary() const;
};

struct DetectorEdgeData {
    DetectorEdgeId detectors;
    std::vector<size_t> observables;     /// The observables crossed along this edge
    double weight;                       /// The weight of the edge to this neighboring node
    double error_probability;            /// The error probability associated with this node
    double correlated_proabilities_sum;  /// The error probability associated with this node
    std::vector<ImpliedWeightUnconverted> implied_weights_for_other_edges{};
    bool operator==(const DetectorEdgeData &other) const;
    bool operator!=(const DetectorEdgeData &other) const;
};

struct DecomposedDemError {
    /// The probability of this error occurring.
    double probability;
    /// Effects of the error. In normal surface code we never decompose into more
    /// than 4, but here it's 8 as a safety factor. It's generally not a great
    /// sign to see so many components, but it's more frustrating to have the
    /// decoding infrastructure refuse to even accept a wonky DEM.
    stim::FixedCapVector<DetectorEdgeData, 8> components;

    bool operator==(const DecomposedDemError &other) const;
    bool operator!=(const DecomposedDemError &other) const;
};

double bernoulli_xor(double p1, double p2);

void add_decomposed_error_to_conditional_groups(
    pm::DecomposedDemError &error, std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> &conditional_groups);

}  // namespace pm
#endif
