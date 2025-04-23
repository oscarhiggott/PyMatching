// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "pymatching/sparse_blossom/driver/io.h"

#include <algorithm>

namespace {

double bernoulli_xor(double p1, double p2) {
    return p1 * (1 - p2) + p2 * (1 - p1);
}

}  // namespace

pm::DetectorEdgeId::DetectorEdgeId() : d1(SIZE_MAX), d2(SIZE_MAX) {
}

pm::DetectorEdgeId::DetectorEdgeId(size_t d) : d1(d), d2(SIZE_MAX) {
}

pm::DetectorEdgeId::DetectorEdgeId(size_t node1, size_t node2) : d1(node1), d2(node2) {
    if (d2 < d1) {
        std::swap(d1, d2);
    }
}

bool pm::DetectorEdgeId::operator==(const DetectorEdgeId &other) const {
    return d1 == other.d1 && d2 == other.d2;
}

bool pm::DetectorEdgeId::operator!=(const DetectorEdgeId &other) const {
    return !(*this == other);
}

bool pm::DetectorEdgeId::operator<(const DetectorEdgeId &other) const {
    return d1 != other.d1 ? d1 < other.d1 : d2 < other.d2;
}

bool pm::DetectorEdgeId::is_boundary() const {
    return d2 == SIZE_MAX;
}

bool pm::DetectorEdgeId::is_valid_edge() const {
    return d1 != SIZE_MAX;
}

bool pm::DetectorEdgeData::operator==(const DetectorEdgeData &other) const {
    return detectors == other.detectors && weight == other.weight && error_probability == other.error_probability &&
           observables == other.observables;
}

bool pm::DetectorEdgeData::operator!=(const DetectorEdgeData &other) const {
    return !(*this == other);
}
bool pm::DecomposedDemError::operator==(const DecomposedDemError &other) const {
    return (components == other.components) && (probability == other.probability);
}

bool pm::DecomposedDemError::operator!=(const DecomposedDemError &other) const {
    return !(*this == other);
}

double pm::merge_weights(double a, double b) {
    auto sgn = std::copysign(1, a) * std::copysign(1, b);
    auto signed_min = sgn * std::min(std::abs(a), std::abs(b));
    return signed_min + std::log(1 + std::exp(-std::abs(a + b))) - std::log(1 + std::exp(-std::abs(a - b)));
}

void pm::add_decomposed_error_to_conditional_groups(
    pm::DecomposedDemError &error, std::map<DetectorEdgeId, std::map<DetectorEdgeId, double>> &conditional_groups) {
    if (error.components.size() == 1) {
        auto c = error.components[0].detectors;
        auto &p = conditional_groups[c][c];
        p = bernoulli_xor(p, error.probability);
    } else {
        for (size_t k1 = 0; k1 < error.components.size(); k1++) {
            for (size_t k2 = k1 + 1; k2 < error.components.size(); k2++) {
                auto c1 = error.components[k1].detectors;
                auto c2 = error.components[k2].detectors;
                auto &p1 = conditional_groups[c1][c2];
                auto &p2 = conditional_groups[c2][c1];
                p1 = bernoulli_xor(p1, error.probability);
                p2 = bernoulli_xor(p2, error.probability);
            }
        }
    }
}
