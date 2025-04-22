#include "pymatching/sparse_blossom/driver/io_helpers.h"

namespace pm {

DetectorEdgeNodes::DetectorEdgeNodes() : d1(SIZE_MAX), d2(SIZE_MAX) {
}

DetectorEdgeNodes::DetectorEdgeNodes(size_t d) : d1(d), d2(SIZE_MAX) {
}

DetectorEdgeNodes::DetectorEdgeNodes(size_t node1, size_t node2) : d1(node1), d2(node2) {
    if (d2 < d1) {
        std::swap(d1, d2);
    }
}

bool DetectorEdgeNodes::operator==(const DetectorEdgeNodes &other) const {
    return d1 == other.d1 && d2 == other.d2;
}

bool DetectorEdgeNodes::operator!=(const DetectorEdgeNodes &other) const {
    return !(*this == other);
}

bool DetectorEdgeNodes::operator<(const DetectorEdgeNodes &other) const {
    return d1 != other.d1 ? d1 < other.d1 : d2 < other.d2;
}

bool DetectorEdgeNodes::is_boundary() const {
    return d2 == SIZE_MAX;
}

bool UserEdge::operator==(const UserEdge &other) const {
    return detectors == other.detectors && weight == other.weight && error_probability == other.error_probability &&
           observables == other.observables;
}

bool UserEdge::operator!=(const UserEdge &other) const {
    return !(*this == other);
}

bool DecomposedDemError::operator==(const DecomposedDemError &other) const {
    return (components == other.components) && (probability == other.probability);
}

bool DecomposedDemError::operator!=(const DecomposedDemError &other) const {
    return !(*this == other);
}

double bernoulli_xor_tmp(double p1, double p2) {
    return p1 * (1 - p2) + p2 * (1 - p1);
}

void add_decomposed_error_to_conditional_groups(
    pm::DecomposedDemError &error,
    std::map<DetectorEdgeNodes, std::map<DetectorEdgeNodes, double>> &conditional_groups) {
    if (error.components.size() == 1) {
        auto c = error.components[0].detectors;
        auto &p = conditional_groups[c][c];
        p = bernoulli_xor_tmp(p, error.probability);
    } else {
        for (size_t k1 = 0; k1 < error.components.size(); k1++) {
            for (size_t k2 = k1 + 1; k2 < error.components.size(); k2++) {
                auto c1 = error.components[k1].detectors;
                auto c2 = error.components[k2].detectors;
                auto &p1 = conditional_groups[c1][c2];
                auto &p2 = conditional_groups[c2][c1];
                p1 = bernoulli_xor_tmp(p1, error.probability);
                p2 = bernoulli_xor_tmp(p2, error.probability);
            }
        }
    }
}

}  // namespace pm
