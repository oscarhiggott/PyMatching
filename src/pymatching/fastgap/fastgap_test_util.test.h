// Copyright 2026 fastgap contributors
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

#ifndef PYMATCHING2_FASTGAP_TEST_UTIL_TEST_H
#define PYMATCHING2_FASTGAP_TEST_UTIL_TEST_H

#include <random>
#include <set>
#include <sstream>
#include <vector>

#include "pymatching/fastgap/gap_index.h"
#include "stim.h"

namespace pm {
namespace fastgap {
namespace testutil {

/// A random matchable DEM: a w x h grid (with optional random diagonals, which create odd cycles
/// and therefore blossoms), boundary edges on the left (observable side) and right columns, random
/// probabilities, random duplicate mechanisms, and observable labels scrambled by a random gauge so
/// that the observable appears on interior edges too.
inline stim::DetectorErrorModel random_grid_dem(
    std::mt19937_64& rng, size_t w, size_t h, double diag_frac = 0.3, double p_lo = 0.01, double p_hi = 0.2) {
    std::uniform_real_distribution<double> up(p_lo, p_hi), u01(0, 1);
    size_t n = w * h;
    std::vector<uint8_t> s(n);
    for (auto& b : s)
        b = rng() & 1;
    auto id = [&](size_t x, size_t y) { return y * w + x; };
    std::ostringstream out;
    auto edge = [&](size_t a, size_t b, int obs) {
        int o = obs ^ s[a] ^ s[b];
        int copies = u01(rng) < 0.15 ? 2 : 1;
        for (int c = 0; c < copies; c++) {
            out << "error(" << up(rng) << ") D" << a << " D" << b;
            if (o)
                out << " L0";
            out << "\n";
        }
    };
    auto bedge = [&](size_t a, int obs) {
        int o = obs ^ s[a];
        int copies = u01(rng) < 0.15 ? 2 : 1;
        for (int c = 0; c < copies; c++) {
            out << "error(" << up(rng) << ") D" << a;
            if (o)
                out << " L0";
            out << "\n";
        }
    };
    for (size_t y = 0; y < h; y++) {
        for (size_t x = 0; x < w; x++) {
            if (x + 1 < w)
                edge(id(x, y), id(x + 1, y), 0);
            if (y + 1 < h)
                edge(id(x, y), id(x, y + 1), 0);
            if (x + 1 < w && y + 1 < h && u01(rng) < diag_frac)
                edge(id(x, y), id(x + 1, y + 1), 0);
        }
        bedge(id(0, y), 1);
        bedge(id(w - 1, y), 0);
    }
    for (size_t i = 0; i < n; i++)
        out << "detector(" << (i % w) << ", " << (i / w) << ") D" << i << "\n";
    return stim::DetectorErrorModel(out.str());
}

/// Exact sample of a DEM's detection events: each error instruction fires independently.
inline std::vector<uint64_t> sample_dem(const stim::DetectorErrorModel& dem, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> u01(0, 1);
    std::set<uint64_t> on;
    dem.iter_flatten_error_instructions([&](const stim::DemInstruction& e) {
        if (u01(rng) < e.arg_data[0]) {
            for (auto& t : e.target_data) {
                if (t.is_relative_detector_id()) {
                    auto [it, inserted] = on.insert(t.val());
                    if (!inserted)
                        on.erase(it);
                }
            }
        }
    });
    return {on.begin(), on.end()};
}

inline stim::DetectorErrorModel surface_code_dem(uint32_t d, uint64_t rounds, double p, const std::string& task) {
    stim::CircuitGenParameters params(rounds, d, task);
    params.after_clifford_depolarization = p;
    params.before_round_data_depolarization = p;
    params.before_measure_flip_probability = p;
    params.after_reset_flip_probability = p;
    auto circuit = stim::generate_surface_code_circuit(params).circuit;
    return stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, true, true, false, 1, false, false);
}

/// Index with a few single-node landmarks spread over the graph (enough to exercise LB).
inline std::shared_ptr<GapIndex> make_index(
    const stim::DetectorErrorModel& dem, size_t num_landmarks = 4) {
    auto idx = GapIndex::create_from_dem(dem.str(), 0);
    std::vector<LandmarkSource> sources;
    size_t n = idx->graph.num_nodes;
    for (size_t i = 0; i < num_landmarks && n > 0; i++)
        sources.push_back({"node" + std::to_string(i), {(uint32_t)((i * 7919) % n)}});
    idx->set_landmarks(sources, "{\"strategy\":\"test\"}");
    return idx;
}

}  // namespace testutil
}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_TEST_UTIL_TEST_H
