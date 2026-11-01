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

#ifndef PYMATCHING2_FASTGAP_GAP_DECODER_H
#define PYMATCHING2_FASTGAP_GAP_DECODER_H

#include <algorithm>
#include <atomic>
#include <memory>
#include <vector>

#include "pymatching/fastgap/dual_state.h"
#include "pymatching/fastgap/gap_index.h"
#include "pymatching/fastgap/thread_pool.h"

namespace pm {
namespace fastgap {

enum class Pricing : uint8_t {
    /// sigma_hat from the landmark bound (CLAUDE.md §3.4): the number fastgap reports.
    INDEX = 0,
    /// Exact sigma from a per-shot Dijkstra: the slow reference Delta_sigma (CLAUDE.md §5.2).
    EXACT = 1,
};

struct DecodeOptions {
    /// Early-termination threshold tau in internal units (CLAUDE.md §3.6); < 0 means none.
    dist_int threshold = -1;
    /// Compute gap_ub (CLAUDE.md §3.5).
    bool upper_bound = false;
    Pricing pricing = Pricing::INDEX;
    /// Run the §5.1 debug invariants on every shot (always on in tests).
    bool check = CHECKS_DEFAULT_ON;
    /// Keep the walk in the result (single-shot API).
    bool keep_walk = false;
};

/// One step of the returned walk P_H: a hop from `from` to `to`. Endpoints are detector ids, or
/// -1 for L and -2 for R.
struct Hop {
    int64_t from;
    int64_t to;
    dist_int price;  // c(from, to) used by the search
    bool exact;      // the price used an exact distance
};

struct GapResult {
    pm::obs_int prediction = 0;  // original (un-gauged) observable labels
    dist_int w_star = 0;
    dist_int gap_lb = 0;
    dist_int gap_ub = DIST_INF;
    bool censored = false;
    bool walk_simple = false;
    bool all_hops_exact = false;
    bool blossom_pass = false;  // some hop-in/hop-out run through a blossom
    uint32_t num_defects = 0;
    uint32_t num_blossoms = 0;
    uint32_t num_hops = 0;
    dist_int hop_slack = -1;  // sum over hops of (dist - d_hat), only with upper_bound
    int64_t t_decode_ns = 0;
    int64_t t_gap_ns = 0;
    std::vector<Hop> walk;
};

/// Decoder state: the Mwpm plus scratch buffers reused across shots.
struct GapWorker {
    pm::Mwpm mwpm;
    DualState ds;
    DualStateScratch ds_scratch;
    std::vector<uint64_t> events;
    std::vector<dist_int> label;
    std::vector<uint8_t> settled;
    std::vector<int32_t> pred_state;
    std::vector<int32_t> pred_hop;
    std::vector<dist_int> row;  // exact pricing: distances from the current state to all detectors
    InteriorGraph::AStarScratch astar;
    // Per-shot pricing scratch.
    std::vector<int32_t> act;             // exact pricing: active hop targets (defect indices), compacted
    std::vector<int32_t> pos;             // exact pricing: position of each defect in `act`, or -1
    std::vector<int32_t> top;             // top-level blossom of each defect, -1 if none
};

/// State of one shot's search with index pricing (CLAUDE.md §3.3, §3.4).
///
/// Every hop target u (a defect not matched to L) leads, by its relay, to exactly one place:
/// the state mate(u), or R. Labels are therefore stored per target position, which makes
/// relaxing a row one contiguous, branch-free pass. Targets leading to states come first,
/// sorted by state index (so the first minimum is the serial tie-break), then targets leading to
/// R, sorted by defect index. The start states (L, and defects matched to L) all have label 0
/// and are never relaxed, so they are kept in a separate sorted list.
///
/// With several threads, thread 0 runs the search while the others prefetch price rows "standing
/// at x" for the unsettled states with the smallest labels, the ones thread 0 settles next. A
/// row's content does not depend on who priced it, so results never depend on the thread count.
struct RowCache {
    static constexpr uint32_t ROW_FREE = 0, ROW_CLAIMED = 1, ROW_READY = 2;
    size_t num_states = 0;  // n + 1 (state n is "standing at L")
    size_t num_targets = 0;
    size_t num_state_targets = 0;   // targets [0, num_state_targets) lead to a state
    std::vector<int32_t> tgt_u;     // target defect, per position
    std::vector<int32_t> tgt_state; // state reached by the relay, per position < num_state_targets
    std::vector<int32_t> tgt_pos;   // position of each defect among the targets, or -1
    std::vector<int32_t> starts;    // start states, ascending
    std::vector<int32_t> gathered;  // target potentials, potential-major (k * num_targets + p)
    // Thread 0's search state, per target position.
    std::vector<dist_int> label;
    std::vector<dist_int> settled_mask;  // 0, or DIST_INF once the state is settled
    std::vector<int32_t> pred_state;
    std::vector<int32_t> pred_hop;
    std::vector<dist_int> hint_shadow;  // last label published to `hint`, per target position
    // Shared with the prefetchers. Rows are stored as int32, saturated at `cap` = d_LR + 1: the
    // answer is at most d_LR, so a state whose label reaches d_LR + 1 is never settled, never on
    // the walk, and never the answer, and saturation changes nothing that is reported.
    dist_int cap = DIST_INF;
    std::vector<int32_t> rows;  // num_states x num_targets prices
    std::unique_ptr<std::atomic<uint32_t>[]> row_state;
    std::unique_ptr<std::atomic<dist_int>[]> hint;  // label published by thread 0 (INF once settled)
    size_t capacity = 0;
    std::atomic<bool> done{false};

    void reserve_states(size_t n) {
        if (n <= capacity)
            return;
        capacity = std::max(n, 2 * capacity);
        row_state = std::make_unique<std::atomic<uint32_t>[]>(capacity);
        hint = std::make_unique<std::atomic<dist_int>[]>(capacity);
    }
};

/// Per-thread scratch for pricing one row.
struct alignas(64) RowScratch {
    std::vector<int32_t> lb;
    std::vector<dist_int> row;  // single-threaded search: the current row
};

class GapDecoder {
   public:
    /// `num_threads` threads cooperate on each shot's gap search (intra-shot parallelism); shots
    /// are always decoded one after another.
    GapDecoder(std::shared_ptr<const GapIndex> index, size_t num_threads = 1, bool pin_threads = false);

    /// Decode one shot given as detection-event indices.
    GapResult decode(const std::vector<uint64_t>& detection_events, const DecodeOptions& options);

    /// Decode `num_shots` bit-packed shots (`bytes_per_shot` bytes each, little-endian bit order as
    /// produced by stim), in order. `out` is resized to num_shots.
    void decode_batch(
        const uint8_t* shots,
        size_t num_shots,
        size_t bytes_per_shot,
        const DecodeOptions& options,
        std::vector<GapResult>& out);

    const GapIndex& index() const {
        return *index_;
    }
    size_t num_threads() const {
        return team_.size();
    }

    /// Shots with fewer than `parallel_grain` defects are searched on one thread (the row
    /// prefetchers cannot get ahead on tiny shots). Results never depend on it.
    size_t parallel_grain = 16;

   private:
    void decode_with(GapWorker& w, const std::vector<uint64_t>& events, const DecodeOptions& o, GapResult& res);
    std::shared_ptr<const GapIndex> index_;
    SpinTeam team_;
    GapWorker worker_;
    RowCache cache_;
    std::vector<std::unique_ptr<RowScratch>> scratch_;
};

/// Unpacks one bit-packed shot into detection events (dropping `ignored_detector` if given).
void unpack_shot(
    const uint8_t* shot, size_t bytes_per_shot, size_t num_detectors, size_t ignored_detector,
    std::vector<uint64_t>& out);

}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_GAP_DECODER_H
