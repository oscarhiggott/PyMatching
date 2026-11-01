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

#include "pymatching/fastgap/gap_decoder.h"

#include <algorithm>
#include <chrono>
#include <limits>
#include <type_traits>

namespace pm {
namespace fastgap {

namespace {

constexpr int32_t STATE_L = -1;
constexpr int32_t END_R = -2;
constexpr int32_t PRED_FREE_RELAY = -3;  // standing at u with mate(u) = L, reached for free

inline int64_t elapsed_ns(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
}

/// Hop prices c(x, u) (CLAUDE.md §3.3, §3.4). Every price is a lower bound on the true reduced
/// cost, clamped at 0 (safe because sigma >= 0). `hop`, `from_l`, `to_r` and `reduce` are
/// read-only and safe to call from several threads.
///
/// Index pricing prices whole rows in `compute_row` (see RowCache). Exact pricing, the slow
/// reference Delta_sigma, reads one Dijkstra row at a time ("standing at x", `begin_row`) over the
/// *active* hop targets: defects u whose mate is R or an unsettled state (once state m is settled,
/// hopping into mate(m) can no longer improve anything).
struct Pricer {
    const GapIndex& idx;
    const DualState& ds;
    GapWorker& w;
    Pricing mode;
    int32_t row_state = INT32_MIN;
    size_t n = 0;

    Pricer(const GapIndex& idx, const DualState& ds, GapWorker& w, Pricing mode)
        : idx(idx), ds(ds), w(w), mode(mode), n(ds.num_defects()) {
        w.top.resize(n);
        for (size_t i = 0; i < n; i++)
            w.top[i] = ds.top_blossom(i);
        if (mode == Pricing::INDEX)
            return;
        // Hops into a defect matched to L only lead back to L, so those are never targets.
        w.act.clear();
        w.pos.assign(n, -1);
        for (size_t i = 0; i < n; i++) {
            if (ds.mate[i] != MATE_L) {
                w.pos[i] = (int32_t)w.act.size();
                w.act.push_back((int32_t)i);
            }
        }
    }
    Pricer(const Pricer&) = delete;

    size_t num_active() const {
        return w.act.size();
    }

    /// Removes u from the active targets (swap with the last one). Exact pricing only.
    void deactivate(int32_t u) {
        int32_t p = w.pos[u];
        if (p < 0)
            return;
        int32_t moved = w.act.back();
        w.act[p] = moved;
        w.pos[moved] = p;
        w.pos[u] = -1;
        w.act.pop_back();
    }

    /// Exact pricing only: distances from state x to every detector.
    void begin_row(int32_t x) {
        if (x < 0 || row_state == x)
            return;
        row_state = x;
        w.row = idx.graph.dijkstra({{ds.defects[x], 0}});
    }

    /// Exact pricing only: price of the hop from the current row x to the active target at
    /// position p (requires begin_row(x)).
    inline dist_int price_at(int32_t x, size_t p) const {
        int32_t j = w.act[p];
        return reduce(x, j, w.row[ds.defects[j]]);
    }

    inline dist_int beta(int32_t i, int32_t j) const {
        return (w.top[i] < 0 || w.top[i] != w.top[j]) ? 0 : ds.beta(i, j);
    }

    inline dist_int reduce(int32_t i, int32_t j, dist_int d) const {
        if (is_inf(d))
            return DIST_INF;
        dist_int s = d - ds.rho[i] - ds.rho[j] + 2 * beta(i, j);
        return s > 0 ? s : 0;
    }

    /// Price of any hop i -> j; read-only, so safe to call from several threads in INDEX mode.
    /// In EXACT mode it requires begin_row(i). Optionally returns the distance estimate.
    inline dist_int hop(int32_t i, int32_t j, bool& exact, dist_int* d_out = nullptr) const {
        dist_int d;
        if (mode == Pricing::EXACT) {
            exact = true;
            d = w.row[ds.defects[j]];
        } else {
            d = idx.d_hat(ds.defects[i], ds.defects[j], exact);
        }
        if (d_out)
            *d_out = d;
        return reduce(i, j, d);
    }

    inline dist_int from_l(int32_t j) const {
        dist_int h = idx.h_l[ds.defects[j]];
        if (is_inf(h))
            return DIST_INF;
        dist_int s = h - ds.rho[j];
        return s > 0 ? s : 0;
    }

    inline dist_int to_r(int32_t i) const {
        if (i == STATE_L)
            return idx.d_lr;
        dist_int h = idx.h_r[ds.defects[i]];
        if (is_inf(h))
            return DIST_INF;
        dist_int s = h - ds.rho[i];
        return s > 0 ? s : 0;
    }
};

struct InternalHop {
    int32_t from;  // defect index or STATE_L
    int32_t to;    // defect index or END_R
    dist_int price;
    dist_int dhat;  // the distance estimate behind the price
    bool exact;
};

/// Debug invariants of CLAUDE.md §5.1: relays tight, sigma >= 0, d_hat <= dist.
void check_dual_invariants(const GapIndex& idx, const DualState& ds, GapWorker& w) {
    size_t n = ds.num_defects();
    for (size_t i = 0; i < n; i++) {
        uint32_t u = ds.defects[i];
        if (idx.h_l[u] - ds.rho[i] < 0 || idx.h_r[u] - ds.rho[i] < 0)
            throw InvariantViolation("boundary reduced cost sigma(u, X) < 0 at D" + std::to_string(u));
        int32_t m = ds.mate[i];
        if (m == MATE_L && idx.h_l[u] != ds.rho[i])
            throw InvariantViolation("relay to L is not tight at D" + std::to_string(u));
        if (m == MATE_R && idx.h_r[u] != ds.rho[i])
            throw InvariantViolation("relay to R is not tight at D" + std::to_string(u));
    }
    size_t pair_rows = std::min<size_t>(n, 64);
    std::vector<uint8_t> relay_checked(n, 0);
    for (size_t i = 0; i < n; i++) {
        int32_t m = ds.mate[i];
        bool need_row = i < pair_rows || (m >= 0 && (size_t)m > i);
        if (!need_row)
            continue;
        auto row = idx.graph.dijkstra({{ds.defects[i], 0}});
        for (size_t j = i + 1; j < n; j++) {
            bool is_mate = (int32_t)j == m;
            if (i >= pair_rows && !is_mate)
                continue;
            dist_int d = row[ds.defects[j]];
            dist_int reduced = d - ds.rho[i] - ds.rho[j] + 2 * ds.beta(i, j);
            if (!is_inf(d) && reduced < 0)
                throw InvariantViolation(
                    "dual infeasible: sigma(D" + std::to_string(ds.defects[i]) + ", D" + std::to_string(ds.defects[j]) +
                    ") = " + std::to_string(reduced) + " < 0");
            if (is_mate) {
                if (reduced != 0)
                    throw InvariantViolation("relay between matched defects is not tight");
                relay_checked[i] = relay_checked[j] = 1;
            }
            bool exact;
            dist_int dh = idx.d_hat(ds.defects[i], ds.defects[j], exact);
            if (dh > d || (exact && dh != d))
                throw InvariantViolation("d_hat overestimates the interior distance");
        }
    }
    for (size_t i = 0; i < n; i++) {
        if (ds.mate[i] >= 0 && !relay_checked[i])
            throw InvariantViolation("relay tightness was not checked");
    }
}

}  // namespace

void unpack_shot(
    const uint8_t* shot, size_t bytes_per_shot, size_t num_detectors, size_t ignored_detector,
    std::vector<uint64_t>& out) {
    out.clear();
    size_t limit = std::min(num_detectors, bytes_per_shot * 8);
    for (size_t b = 0; b < bytes_per_shot; b++) {
        uint8_t byte = shot[b];
        while (byte) {
            int bit = __builtin_ctz(byte);
            byte &= byte - 1;
            size_t d = b * 8 + bit;
            if (d >= limit)
                throw std::invalid_argument(
                    "shot has a detection event at index " + std::to_string(d) + " but the graph only has " +
                    std::to_string(num_detectors) + " detectors");
            if (d != ignored_detector)
                out.push_back(d);
        }
    }
}

namespace {

/// Result of the search over "standing at x" states. `best_hop` is END_R (hop straight to R) or
/// the defect whose mate is R; ties are broken by (iteration, hop) so the walk does not depend on
/// the number of threads.
struct SearchOutcome {
    dist_int best = DIST_INF;
    int32_t best_from = INT32_MIN;
    int32_t best_hop = END_R;
    dist_int stop_label = DIST_INF;
};

inline bool better(dist_int c, int64_t it, int32_t hop, dist_int best, int64_t best_it, int32_t best_hop) {
    if (c != best)
        return c < best;
    return it < best_it || (it == best_it && hop < best_hop);
}

/// Dense Dijkstra with exact pricing (CLAUDE.md §3.3, §5.2). State index n is "standing at L".
SearchOutcome exact_search(const DualState& ds, GapWorker& w, Pricer& pricer, dist_int tau) {
    const int32_t n = (int32_t)ds.num_defects();
    SearchOutcome out;
    int64_t best_iter = -1;
    for (int64_t iter = 0;; iter++) {
        int32_t x = -1;
        dist_int lx = DIST_INF;
        for (int32_t s = 0; s <= n; s++) {
            if (!w.settled[s] && w.label[s] < lx) {
                lx = w.label[s];
                x = s;
            }
        }
        if (x < 0 || lx >= out.best || lx >= tau) {
            out.stop_label = lx;
            return out;
        }
        w.settled[x] = 1;
        int32_t xs = x == n ? STATE_L : x;
        // Move B: hop to R.
        dist_int c_r = pricer.to_r(xs);
        if (!is_inf(c_r) && better(lx + c_r, iter, END_R, out.best, best_iter, out.best_hop)) {
            out.best = lx + c_r;
            out.best_from = xs;
            out.best_hop = END_R;
            best_iter = iter;
        }
        // Move A: hop to defect u, then relay to mate(u). Settling x retires the target mate(x).
        if (xs >= 0 && ds.mate[xs] >= 0)
            pricer.deactivate(ds.mate[xs]);
        pricer.begin_row(xs);
        size_t a = pricer.num_active();
        for (size_t p = 0; p < a; p++) {
            int32_t u = w.act[p];
            if (u == xs)
                continue;
            dist_int c = xs == STATE_L ? pricer.from_l(u) : pricer.price_at(xs, p);
            if (is_inf(c))
                continue;
            dist_int cand = lx + c;
            int32_t m = ds.mate[u];
            if (m == MATE_R) {
                if (better(cand, iter, u, out.best, best_iter, out.best_hop)) {
                    out.best = cand;
                    out.best_from = xs;
                    out.best_hop = u;
                    best_iter = iter;
                }
            } else if (cand < w.label[m]) {
                w.label[m] = cand;
                w.pred_state[m] = xs;
                w.pred_hop[m] = u;
            }
        }
    }
}

/// Lays out the targets of this shot (see RowCache) and resets the search state.
void setup_row_cache(const GapIndex& idx, const DualState& ds, RowCache& c, bool shared) {
    const int32_t n = (int32_t)ds.num_defects();
    const size_t stride = idx.potential_stride();
    c.num_states = (size_t)n + 1;
    c.tgt_u.clear();
    c.tgt_state.clear();
    c.starts.clear();
    // Targets leading to a state, by state index; then targets leading to R, by defect index.
    for (int32_t v = 0; v < n; v++) {
        int32_t m = ds.mate[v];
        if (m == MATE_L)
            c.starts.push_back(v);
        else if (m >= 0) {
            // v is the state reached from target m.
            c.tgt_u.push_back(m);
            c.tgt_state.push_back(v);
        }
    }
    c.starts.push_back(n);
    c.num_state_targets = c.tgt_u.size();
    for (int32_t u = 0; u < n; u++) {
        if (ds.mate[u] == MATE_R)
            c.tgt_u.push_back(u);
    }
    const size_t a = c.tgt_u.size();
    const size_t as = c.num_state_targets;
    c.num_targets = a;
    c.tgt_pos.assign(n, -1);
    for (size_t p = 0; p < a; p++)
        c.tgt_pos[c.tgt_u[p]] = (int32_t)p;
    c.gathered.resize(a * stride);
    for (size_t p = 0; p < a; p++) {
        const int32_t* src = idx.potentials32.data() + (size_t)ds.defects[c.tgt_u[p]] * stride;
        for (size_t k = 0; k < stride; k++)
            c.gathered[k * a + p] = src[k];
    }
    c.label.assign(as, DIST_INF);
    c.settled_mask.assign(as, 0);
    c.pred_state.assign(as, INT32_MIN);
    c.pred_hop.assign(as, -1);
    if (!shared)
        return;
    c.cap = idx.d_lr + 1;
    c.hint_shadow.assign(as, DIST_INF);
    c.rows.resize(c.num_states * a);
    c.reserve_states(c.num_states);
    for (size_t s = 0; s < c.num_states; s++) {
        c.row_state[s].store(RowCache::ROW_FREE, std::memory_order_relaxed);
        c.hint[s].store(DIST_INF, std::memory_order_relaxed);
    }
    for (int32_t s : c.starts)
        c.hint[s].store(0, std::memory_order_relaxed);
    c.done.store(false, std::memory_order_relaxed);
}

/// Prices the hops from state x (n is L) to every target into `row`, saturated at `cap`; the
/// hop from x to itself is not allowed and gets `cap`. Same prices as Pricer::begin_row.
template <typename RowT>
void compute_row(const GapIndex& idx, const DualState& ds, const Pricer& pricer, const GapWorker& w,
                 const RowCache& c, RowScratch& sc, int32_t x, RowT* row, dist_int cap) {
    const int32_t n = (int32_t)ds.num_defects();
    const size_t a = c.num_targets;
    if (x == n) {
        for (size_t p = 0; p < a; p++)
            row[p] = (RowT)std::min(pricer.from_l(c.tgt_u[p]), cap);
        return;
    }
    const size_t stride = idx.potential_stride();
    sc.lb.resize(a);
    int32_t* lb32 = sc.lb.data();
    std::fill(lb32, lb32 + a, 0);
    uint32_t ux = ds.defects[x];
    const int32_t* px = idx.potentials32.data() + (size_t)ux * stride;
    for (size_t k = 0; k < stride; k++) {
        const int32_t* col = c.gathered.data() + k * a;
        int32_t pk = px[k];
        for (size_t p = 0; p < a; p++) {
            int32_t diff = col[p] - pk;
            diff = diff < 0 ? -diff : diff;
            lb32[p] = diff > lb32[p] ? diff : lb32[p];
        }
    }
    for (size_t p = 0; p < a; p++)
        row[p] = (RowT)std::min(pricer.reduce(x, c.tgt_u[p], lb32[p]), cap);
    if (c.tgt_pos[x] >= 0)
        row[c.tgt_pos[x]] = (RowT)cap;
}

/// Prefetching thread: repeatedly claims the free state with the smallest published label and
/// prices its row, until thread 0 finishes the search.
void prefetch_rows(const GapIndex& idx, const DualState& ds, const Pricer& pricer, const GapWorker& w, RowCache& c,
                   RowScratch& sc) {
    const size_t ns = c.num_states;
    const size_t a = c.num_targets;
    while (!c.done.load(std::memory_order_acquire)) {
        int32_t pick = -1;
        dist_int best = DIST_INF;
        for (size_t s = 0; s < ns; s++) {
            dist_int l = c.hint[s].load(std::memory_order_relaxed);
            if (l < best && c.row_state[s].load(std::memory_order_relaxed) == RowCache::ROW_FREE) {
                best = l;
                pick = (int32_t)s;
            }
        }
        if (pick < 0) {
            cpu_relax();
            continue;
        }
        uint32_t expected = RowCache::ROW_FREE;
        if (!c.row_state[pick].compare_exchange_strong(expected, RowCache::ROW_CLAIMED, std::memory_order_acq_rel))
            continue;
        compute_row(idx, ds, pricer, w, c, sc, pick, c.rows.data() + (size_t)pick * a, c.cap);
        c.row_state[pick].store(RowCache::ROW_READY, std::memory_order_release);
    }
}

/// Dense Dijkstra with index pricing: the same iterations, prices and tie-breaks as
/// exact_search (exact pricing). With `shared`, rows come from the prefetchers when ready.
template <bool shared>
SearchOutcome row_search(const GapIndex& idx, const DualState& ds, const Pricer& pricer, GapWorker& w, RowCache& c,
                         RowScratch& sc, dist_int tau) {
    using RowT = typename std::conditional<shared, int32_t, dist_int>::type;
    const int32_t n = (int32_t)ds.num_defects();
    const size_t a = c.num_targets;
    const size_t as = c.num_state_targets;
    dist_int* label = c.label.data();
    dist_int* mask = c.settled_mask.data();
    int32_t* pst = c.pred_state.data();
    int32_t* phop = c.pred_hop.data();
    const int32_t* tu = c.tgt_u.data();
    if constexpr (!shared)
        sc.row.resize(a);
    // First position holding the smallest unsettled label `m` (the smallest state index).
    auto first_min = [&](dist_int m) {
        size_t p = 0;
        while (std::max(label[p], mask[p]) != m)
            p++;
        return (int32_t)p;
    };
    size_t next_start = 0;
    SearchOutcome out;
    int64_t best_iter = -1;
    // Smallest unsettled label; each relaxation pass computes it for the next iteration.
    dist_int lmin = DIST_INF;
    for (size_t p = 0; p < as; p++)
        lmin = std::min(lmin, std::max(label[p], mask[p]));
    int32_t next_px = -1;
    for (int64_t iter = 0;; iter++) {
        dist_int lx = lmin;
        int32_t px = -1, x = -1;
        if (lx < DIST_INF) {
            px = next_px >= 0 ? next_px : first_min(lx);
            x = c.tgt_state[px];
        }
        if (next_start < c.starts.size() && (lx > 0 || c.starts[next_start] < x)) {
            px = -1;
            x = c.starts[next_start];
            lx = 0;
        }
        if (x < 0 || lx >= out.best || lx >= tau) {
            out.stop_label = x < 0 ? DIST_INF : lx;
            break;
        }
        if (px >= 0)
            mask[px] = DIST_INF;
        else
            next_start++;
        int32_t xs = x == n ? STATE_L : x;
        // Move B: hop to R.
        dist_int c_r = pricer.to_r(xs);
        if (!is_inf(c_r) && better(lx + c_r, iter, END_R, out.best, best_iter, out.best_hop)) {
            out.best = lx + c_r;
            out.best_from = xs;
            out.best_hop = END_R;
            best_iter = iter;
        }
        // Row x: from a prefetcher, priced here, or (rarely) waited for.
        const RowT* row;
        if constexpr (!shared) {
            compute_row(idx, ds, pricer, w, c, sc, x, sc.row.data(), DIST_INF);
            row = sc.row.data();
        } else {
            c.hint[x].store(DIST_INF, std::memory_order_relaxed);
            int32_t* slot = c.rows.data() + (size_t)x * a;
            uint32_t st = c.row_state[x].load(std::memory_order_acquire);
            if (st == RowCache::ROW_FREE &&
                c.row_state[x].compare_exchange_strong(st, RowCache::ROW_CLAIMED, std::memory_order_acq_rel)) {
                compute_row(idx, ds, pricer, w, c, sc, x, slot, c.cap);
                c.row_state[x].store(RowCache::ROW_READY, std::memory_order_relaxed);
            } else {
                while (c.row_state[x].load(std::memory_order_acquire) != RowCache::ROW_READY)
                    cpu_relax();
            }
            row = slot;
        }
        // Move A into a state: hop to u, relay to mate(u). Branch-free; a settled state (mate(x)
        // in particular) never improves, since cand >= lx >= its label; INF prices never do.
        dist_int nm = DIST_INF;
        for (size_t p = 0; p < as; p++) {
            dist_int cand = lx + (dist_int)row[p];
            bool imp = cand < label[p];
            label[p] = imp ? cand : label[p];
            pst[p] = imp ? xs : pst[p];
            phop[p] = imp ? tu[p] : phop[p];
            nm = std::min(nm, std::max(label[p], mask[p]));
        }
        lmin = nm;
        next_px = -1;
        if constexpr (shared) {
            // Pull the next state's row towards this core while the rest of the iteration runs.
            if (nm < DIST_INF) {
                next_px = first_min(nm);
                const int32_t* next_row = c.rows.data() + (size_t)c.tgt_state[next_px] * a;
                for (size_t off = 0; off < a; off += 64 / sizeof(int32_t))
                    __builtin_prefetch(next_row + off);
            }
        }
        // Move A into R: the cheapest target wins, ties to the smallest defect index.
        for (size_t p = as; p < a; p++) {
            if (is_inf((dist_int)row[p]))
                continue;
            dist_int cand = lx + (dist_int)row[p];
            if (better(cand, iter, tu[p], out.best, best_iter, out.best_hop)) {
                out.best = cand;
                out.best_from = xs;
                out.best_hop = tu[p];
                best_iter = iter;
            }
        }
        if constexpr (shared) {
            // Publish the labels that dropped, as the prefetchers' priorities.
            dist_int* shadow = c.hint_shadow.data();
            for (size_t p = 0; p < as; p++) {
                if (label[p] < shadow[p] && !mask[p]) {
                    shadow[p] = label[p];
                    c.hint[c.tgt_state[p]].store(label[p], std::memory_order_relaxed);
                }
            }
        }
    }
    for (size_t p = 0; p < as; p++) {
        if (pst[p] != INT32_MIN) {
            w.pred_state[c.tgt_state[p]] = pst[p];
            w.pred_hop[c.tgt_state[p]] = phop[p];
        }
    }
    return out;
}

}  // namespace

GapDecoder::GapDecoder(std::shared_ptr<const GapIndex> index, size_t num_threads, bool pin_threads)
    : index_(std::move(index)), team_(num_threads, pin_threads) {
    if (!index_->finalized)
        throw std::invalid_argument("GapIndex has no landmarks yet (call set_landmarks first)");
    worker_.mwpm = index_->make_mwpm();
    for (size_t t = 0; t < team_.size(); t++)
        scratch_.push_back(std::make_unique<RowScratch>());
}

GapResult GapDecoder::decode(const std::vector<uint64_t>& detection_events, const DecodeOptions& options) {
    GapResult res;
    auto& w = worker_;
    w.events.clear();
    for (auto d : detection_events) {
        if (d != index_->spec.ref_obs_detector)
            w.events.push_back(d);
    }
    decode_with(w, w.events, options, res);
    return res;
}

void GapDecoder::decode_batch(
    const uint8_t* shots,
    size_t num_shots,
    size_t bytes_per_shot,
    const DecodeOptions& options,
    std::vector<GapResult>& out) {
    out.resize(num_shots);
    size_t num_detectors = index_->num_detectors();
    size_t ignored = index_->spec.ref_obs_detector;
    auto& w = worker_;
    for (size_t shot = 0; shot < num_shots; shot++) {
        unpack_shot(shots + shot * bytes_per_shot, bytes_per_shot, num_detectors, ignored, w.events);
        out[shot] = GapResult();
        decode_with(w, w.events, options, out[shot]);
    }
}

void GapDecoder::decode_with(GapWorker& w, const std::vector<uint64_t>& events, const DecodeOptions& o, GapResult& res) {
    const GapIndex& idx = *index_;
    auto t0 = std::chrono::steady_clock::now();

    // Stage 1: the ordinary, untruncated PyMatching decode, keeping the dual state.
    extract_dual_state(w.mwpm, events, idx.spec.observable, w.ds, w.ds_scratch);
    const DualState& ds = w.ds;
    auto t1 = std::chrono::steady_clock::now();

    pm::obs_int shift = (pm::obs_int)gauge_shift(idx.spec, events) << idx.spec.observable;
    res.prediction = ds.obs_mask ^ shift;
    res.w_star = ds.w_star;
    res.num_defects = (uint32_t)ds.num_defects();
    res.num_blossoms = (uint32_t)ds.num_blossoms();
    if (o.check)
        check_dual_invariants(idx, ds, w);

    // Stage 2: dense Dijkstra over the |D| + 2 states "standing at x" (CLAUDE.md §3.3). With
    // several threads, thread 0 runs it while the others prefetch its price rows.
    const int32_t n = (int32_t)ds.num_defects();
    const bool index_pricing = o.pricing == Pricing::INDEX;
    // The shared rows hold int32 prices saturated at d_LR + 1, which must fit.
    const bool parallel = index_pricing && team_.size() > 1 && (size_t)n >= parallel_grain &&
                          idx.d_lr < (dist_int)std::numeric_limits<int32_t>::max() - 1;
    Pricer pricer(idx, ds, w, o.pricing);

    w.pred_state.assign(n + 1, INT32_MIN);
    w.pred_hop.assign(n + 1, -1);
    for (int32_t j = 0; j < n; j++) {
        if (ds.mate[j] == MATE_L)
            w.pred_state[j] = PRED_FREE_RELAY;
    }
    dist_int tau = o.threshold >= 0 ? o.threshold : DIST_INF;
    SearchOutcome so;
    if (!index_pricing) {
        // State index n is "standing at L".
        w.label.assign(n + 1, DIST_INF);
        w.settled.assign(n + 1, 0);
        w.label[n] = 0;
        for (int32_t j = 0; j < n; j++) {
            if (ds.mate[j] == MATE_L)
                w.label[j] = 0;
        }
        so = exact_search(ds, w, pricer, tau);
    } else if (!parallel) {
        setup_row_cache(idx, ds, cache_, false);
        so = row_search<false>(idx, ds, pricer, w, cache_, *scratch_[0], tau);
    } else {
        setup_row_cache(idx, ds, cache_, true);
        team_.run(team_.size(), [&](size_t t) {
            if (t == 0) {
                so = row_search<true>(idx, ds, pricer, w, cache_, *scratch_[0], tau);
                cache_.done.store(true, std::memory_order_release);
            } else {
                prefetch_rows(idx, ds, pricer, w, cache_, *scratch_[t]);
            }
        });
    }
    dist_int best = so.best;
    int32_t best_from = so.best_from;
    int32_t best_hop = so.best_hop;
    dist_int stop_label = so.stop_label;

    if (best <= stop_label) {
        res.gap_lb = best;
        res.censored = false;
    } else {
        // Every unsettled label is >= tau, so the answer is >= tau (CLAUDE.md §3.6).
        res.gap_lb = tau;
        res.censored = true;
    }

    // Reconstruct the walk P_H (only when it finished).
    std::vector<InternalHop> hops;
    int32_t start_state = STATE_L;
    if (!res.censored && best_from != INT32_MIN) {
        auto make_hop = [&](int32_t a, int32_t b) {
            InternalHop h{a, b, 0, 0, true};
            if (o.pricing == Pricing::EXACT)
                pricer.begin_row(a);
            if (a == STATE_L && b == END_R) {
                h.price = h.dhat = idx.d_lr;
            } else if (a == STATE_L) {
                h.price = pricer.from_l(b);
                h.dhat = idx.h_l[ds.defects[b]];
            } else if (b == END_R) {
                h.price = pricer.to_r(a);
                h.dhat = idx.h_r[ds.defects[a]];
            } else {
                h.price = pricer.hop(a, b, h.exact, &h.dhat);
            }
            return h;
        };
        hops.push_back(make_hop(best_from, best_hop));
        int32_t s = best_from;
        while (s != STATE_L) {
            int32_t p = w.pred_state[s];
            if (p == PRED_FREE_RELAY) {
                start_state = s;
                break;
            }
            hops.push_back(make_hop(p, w.pred_hop[s]));
            s = p;
        }
        std::reverse(hops.begin(), hops.end());
    }

    if (!hops.empty()) {
        // Node sequence z_0 = L, ..., z_last = R, with the type of the edge entering each node.
        struct Step {
            int32_t node;  // defect index, STATE_L or END_R
            bool via_hop;
        };
        std::vector<Step> seq;
        seq.push_back({STATE_L, false});
        if (start_state != STATE_L)
            seq.push_back({start_state, false});
        for (auto& h : hops) {
            seq.push_back({h.to, true});
            if (h.to != END_R) {
                int32_t m = ds.mate[h.to];
                seq.push_back({m >= 0 ? m : END_R, false});
            }
        }
        res.num_hops = (uint32_t)hops.size();
        res.all_hops_exact = std::all_of(hops.begin(), hops.end(), [](const InternalHop& h) { return h.exact; });

        // Simple: every defect at most once.
        std::vector<int32_t> visited;
        for (auto& st : seq) {
            if (st.node >= 0)
                visited.push_back(st.node);
        }
        std::sort(visited.begin(), visited.end());
        res.walk_simple = std::adjacent_find(visited.begin(), visited.end()) == visited.end();

        // Blossom passes: runs inside a blossom S entered and left by hops (n_S of §3.2).
        if (ds.num_blossoms() > 0) {
            std::vector<int32_t> blossoms;
            for (auto v : visited) {
                for (uint32_t a = ds.anc_offset[v]; a < ds.anc_offset[v + 1]; a++)
                    blossoms.push_back(ds.anc_blossom[a]);
            }
            std::sort(blossoms.begin(), blossoms.end());
            blossoms.erase(std::unique(blossoms.begin(), blossoms.end()), blossoms.end());
            auto in_blossom = [&](int32_t v, int32_t b) {
                if (v < 0)
                    return false;
                for (uint32_t a = ds.anc_offset[v]; a < ds.anc_offset[v + 1]; a++) {
                    if (ds.anc_blossom[a] == b)
                        return true;
                }
                return false;
            };
            for (auto b : blossoms) {
                for (size_t p = 0; p < seq.size(); p++) {
                    if (!in_blossom(seq[p].node, b) || (p > 0 && in_blossom(seq[p - 1].node, b)))
                        continue;
                    size_t q = p;
                    while (q + 1 < seq.size() && in_blossom(seq[q + 1].node, b))
                        q++;
                    if (seq[p].via_hop && q + 1 < seq.size() && seq[q + 1].via_hop)
                        res.blossom_pass = true;
                }
            }
        }

        // Upper bound Delta_h on the same walk (CLAUDE.md §3.5).
        if (o.upper_bound) {
            dist_int ub = 0, slack = 0;
            for (auto& h : hops) {
                dist_int d;
                if (h.from == STATE_L && h.to == END_R) {
                    d = idx.d_lr;
                } else if (h.from == STATE_L) {
                    d = idx.h_l[ds.defects[h.to]];
                } else if (h.to == END_R) {
                    d = idx.h_r[ds.defects[h.from]];
                } else if (h.exact) {
                    d = h.dhat;
                } else {
                    d = idx.exact_distance(ds.defects[h.from], ds.defects[h.to], w.astar);
                }
                if (is_inf(d)) {
                    ub = DIST_INF;
                    break;
                }
                if (o.check && d < h.dhat)
                    throw InvariantViolation("d_hat exceeds the exact hop distance");
                slack += d - h.dhat;
                dist_int ra = h.from >= 0 ? ds.r(h.from) : 0;
                dist_int rb = h.to >= 0 ? ds.r(h.to) : 0;
                ub += d - ra - rb;
            }
            res.hop_slack = is_inf(ub) ? -1 : slack;
            res.gap_ub = (res.walk_simple && !is_inf(ub)) ? ub : DIST_INF;
            if (o.check && res.gap_ub < res.gap_lb)
                throw InvariantViolation("gap_ub < gap_lb");
        }

        if (o.keep_walk) {
            auto code = [&](int32_t v) -> int64_t {
                if (v == STATE_L)
                    return -1;
                if (v == END_R)
                    return -2;
                return ds.defects[v];
            };
            for (auto& h : hops)
                res.walk.push_back({code(h.from), code(h.to), h.price, h.exact});
        }
    }
    auto t2 = std::chrono::steady_clock::now();
    res.t_decode_ns = elapsed_ns(t0, t1);
    res.t_gap_ns = elapsed_ns(t1, t2);
}

}  // namespace fastgap
}  // namespace pm
