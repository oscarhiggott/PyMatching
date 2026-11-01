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

#ifndef PYMATCHING2_FASTGAP_GAP_INDEX_H
#define PYMATCHING2_FASTGAP_GAP_INDEX_H

#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "pymatching/fastgap/canonical.h"
#include "pymatching/fastgap/interior_graph.h"

namespace pm {
namespace fastgap {

/// Potentials used for pricing are saturated at POTENTIAL_CAP and stored as int32 (CLAUDE.md
/// §2.1): |min(a, C) - min(b, C)| <= |a - b|, so the bound stays valid, and unreachable entries
/// all become C, so a landmark that reaches neither node contributes 0.
constexpr int32_t POTENTIAL_CAP = std::numeric_limits<int32_t>::max();

inline int32_t saturate_potential(dist_int v) {
    return v >= (dist_int)POTENTIAL_CAP ? POTENTIAL_CAP : (int32_t)v;
}

/// max_i |pu[i] - pv[i]| over two saturated potential rows (branch-free, vectorisable).
inline dist_int lower_bound_rows(const int32_t* pu, const int32_t* pv, size_t k) {
    int32_t best = 0;
    for (size_t i = 0; i < k; i++) {
        int32_t diff = pu[i] - pv[i];  // both in [0, C], so no overflow
        diff = diff < 0 ? -diff : diff;
        best = diff > best ? diff : best;
    }
    return best;
}

/// A landmark: a set of detectors that all sit at distance 0 (CLAUDE.md §7.1).
struct LandmarkSource {
    std::string label;
    std::vector<uint32_t> nodes;
};

/// Per-DEM precomputation (CLAUDE.md §5.3). Built in two steps: `create` canonicalises and builds
/// the graph and boundary potentials; `set_landmarks` adds landmark potentials. Immutable
/// afterwards and shared between decoders and threads.
class GapIndex {
   public:
    GraphSpec spec;
    /// FROM_DEM: the DEM text the index was built from (re-canonicalised on load).
    std::string dem_text;
    InteriorGraph graph;
    /// h_L, h_R: exact distances to the two boundary sides.
    std::vector<dist_int> h_l;
    std::vector<dist_int> h_r;
    /// d_LR = dist(L, R), the weight of the lightest logical.
    dist_int d_lr = DIST_INF;

    /// Landmarks chosen by the strategy (L and R are not in this list; they are always added).
    std::vector<LandmarkSource> landmarks;
    /// Node-major potentials with stride `potential_stride()`: entry 0 is h_L, entry 1 is h_R,
    /// entries 2.. are the landmarks. Exact (int64); `potentials32` is the saturated copy used
    /// for pricing.
    std::vector<dist_int> potentials;
    std::vector<int32_t> potentials32;

    std::string strategy_config;
    uint64_t strategy_hash = 0;
    bool finalized = false;

    static std::shared_ptr<GapIndex> create_from_dem(const std::string& dem_text, size_t observable);
    static std::shared_ptr<GapIndex> create_from_reference_edges(
        size_t num_detectors, size_t obs_detector, std::vector<RefEdge> edges);

    /// Computes the landmark potentials.
    void set_landmarks(std::vector<LandmarkSource> sources, const std::string& config);

    size_t potential_stride() const {
        return landmarks.size() + 2;
    }
    size_t num_detectors() const {
        return spec.num_detectors;
    }

    /// LB(u, v) = max over all (saturated) potentials (h_L, h_R and landmarks) of |h(u) - h(v)|.
    inline dist_int lower_bound(uint32_t u, uint32_t v) const {
        size_t k = potential_stride();
        const int32_t* pu = potentials32.data() + (size_t)u * k;
        const int32_t* pv = potentials32.data() + (size_t)v * k;
        return lower_bound_rows(pu, pv, k);
    }

    /// d_hat(u, v) = LB(u, v) <= dist(u, v). `exact` is set only for u == v; the landmark bound is
    /// never assumed to be exact.
    inline dist_int d_hat(uint32_t u, uint32_t v, bool& exact) const {
        exact = u == v;
        return exact ? 0 : lower_bound(u, v);
    }

    /// Exact interior distance (A* with the landmark bound as heuristic).
    dist_int exact_distance(uint32_t u, uint32_t v, InteriorGraph::AStarScratch& scratch) const;

    /// Creates a fresh Mwpm for decoding (one per thread). Identical to what upstream builds.
    pm::Mwpm make_mwpm() const;

    size_t potential_bytes() const;

    void save(const std::string& path) const;
    static std::shared_ptr<GapIndex> load(const std::string& path);

    /// Recomputes `potentials32` from `potentials`.
    void derive_saturated_potentials();

   private:
    void build_graph();
    mutable std::mutex user_graph_mutex;
    mutable pm::UserGraph decode_user_graph;
};

uint64_t fnv1a(const void* data, size_t len, uint64_t h = 1469598103934665603ULL);

}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_GAP_INDEX_H
