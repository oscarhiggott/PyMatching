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

#include <gtest/gtest.h>

#include <algorithm>

#include "pymatching/fastgap/exact.h"
#include "pymatching/fastgap/fastgap_test_util.test.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"

using namespace pm::fastgap;

namespace {

std::vector<uint8_t> pack(const std::vector<std::vector<uint64_t>>& shots, size_t bytes) {
    std::vector<uint8_t> out(shots.size() * bytes, 0);
    for (size_t s = 0; s < shots.size(); s++) {
        for (auto d : shots[s])
            out[s * bytes + d / 8] |= (uint8_t)(1 << (d % 8));
    }
    return out;
}

struct Counts {
    size_t shots = 0, simple = 0, lb_equals_exact = 0, blossom = 0;
};

/// The full §9.1 chain on one DEM: gap_lb <= Delta_sigma <= Delta == brute force == baseline, and
/// gap_ub >= Delta whenever the walk is simple.
Counts tiny_exact_chain(const stim::DetectorErrorModel& dem, size_t shots, uint64_t seed, size_t num_landmarks) {
    auto idx = testutil::make_index(dem, num_landmarks);
    GapDecoder dec(idx, 1);
    ExactBaseline baseline(idx);
    pm::Mwpm upstream = pm::detector_error_model_to_mwpm(dem, pm::NUM_DISTINCT_WEIGHTS);
    std::mt19937_64 rng(seed);
    DecodeOptions lb_opts;
    lb_opts.upper_bound = true;
    lb_opts.check = true;
    DecodeOptions sigma_opts = lb_opts;
    sigma_opts.pricing = Pricing::EXACT;
    Counts c;
    for (size_t s = 0; s < shots; s++) {
        auto events = testutil::sample_dem(dem, rng);
        if (events.size() > 16)
            continue;
        c.shots++;
        auto lb = dec.decode(events, lb_opts);
        auto sigma = dec.decode(events, sigma_opts);
        auto brute = brute_force_gap(*idx, events);
        auto exact = baseline.decode(events);
        auto ref = pm::decode_detection_events_for_up_to_64_observables(upstream, events, false);

        // Upstream parity (rule 3).
        EXPECT_EQ(lb.w_star, ref.weight);
        EXPECT_EQ(lb.prediction, ref.obs_mask);
        // The exact references agree with each other and with the decoder.
        EXPECT_EQ(exact.gap, brute.gap);
        EXPECT_EQ(exact.w_star, brute.w_star);
        EXPECT_EQ(exact.w_star, lb.w_star);
        if (exact.gap > 0)
            EXPECT_EQ(exact.prediction, lb.prediction & 1);
        // Soundness (rule 1) in integer internal units.
        EXPECT_LE(lb.gap_lb, sigma.gap_lb);
        EXPECT_LE(sigma.gap_lb, exact.gap);
        if (lb.walk_simple) {
            EXPECT_GE(lb.gap_ub, exact.gap);
            c.simple++;
        }
        if (sigma.walk_simple)
            EXPECT_GE(sigma.gap_ub, exact.gap);
        c.lb_equals_exact += lb.gap_lb == exact.gap;
        c.blossom += lb.num_blossoms > 0;
    }
    return c;
}

}  // namespace

TEST(FastgapDecoder, NoDefectsGivesLightestLogical) {
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0 L0
        error(0.1) D0 D1
        error(0.05) D1 D2
        error(0.1) D2
    )DEM");
    auto idx = testutil::make_index(dem);
    GapDecoder dec(idx, 1);
    DecodeOptions o;
    o.upper_bound = true;
    auto r = dec.decode({}, o);
    EXPECT_EQ(r.gap_lb, idx->d_lr);
    EXPECT_EQ(r.gap_ub, idx->d_lr);
    EXPECT_TRUE(r.walk_simple);
    EXPECT_EQ(r.num_hops, 1u);
}

TEST(FastgapDecoder, RepetitionCodeHandComputed) {
    // Chain L - D0 - D1 - D2 - R with equal weights w. One defect at D0: matched to L (weight w);
    // the complement matches it to R (weight 3w), so the gap is 2w.
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0 L0
        error(0.1) D0 D1
        error(0.1) D1 D2
        error(0.1) D2
    )DEM");
    auto idx = testutil::make_index(dem);
    GapDecoder dec(idx, 1);
    DecodeOptions o;
    o.upper_bound = true;
    auto r = dec.decode({0}, o);
    dist_int w = idx->graph.boundary_weight[0];
    EXPECT_EQ(r.w_star, w);
    EXPECT_EQ(r.gap_lb, 2 * w);
    EXPECT_EQ(r.gap_ub, 2 * w);
    EXPECT_EQ(r.prediction, 1u);
}

TEST(FastgapDecoder, TinyExactRandomGrids) {
    std::mt19937_64 rng(1234);
    Counts total;
    for (int trial = 0; trial < 16; trial++) {
        auto dem = testutil::random_grid_dem(rng, 3 + trial % 4, 3 + trial % 3, 0.5, 0.02, 0.25);
        for (size_t k : {(size_t)0, (size_t)3}) {
            auto c = tiny_exact_chain(dem, 80, trial * 31 + k, k);
            total.shots += c.shots;
            total.simple += c.simple;
            total.lb_equals_exact += c.lb_equals_exact;
            total.blossom += c.blossom;
        }
    }
    EXPECT_GT(total.shots, 1000u);
    EXPECT_GT(total.blossom, 0u);
    EXPECT_GT(total.simple, total.shots / 2);
}

TEST(FastgapDecoder, EveryNodeALandmarkGivesExactSigma) {
    std::mt19937_64 rng(77);
    for (int trial = 0; trial < 6; trial++) {
        auto dem = testutil::random_grid_dem(rng, 5, 4, 0.5);
        auto idx = GapIndex::create_from_dem(dem.str(), 0);
        std::vector<LandmarkSource> sources;
        for (uint32_t v = 0; v < idx->graph.num_nodes; v++)
            sources.push_back({"node" + std::to_string(v), {v}});
        idx->set_landmarks(sources, "{}");
        GapDecoder dec(idx, 1);
        DecodeOptions a, b;
        b.pricing = Pricing::EXACT;
        for (int s = 0; s < 100; s++) {
            auto events = testutil::sample_dem(dem, rng);
            auto ra = dec.decode(events, a);
            auto rb = dec.decode(events, b);
            EXPECT_EQ(ra.gap_lb, rb.gap_lb);
        }
    }
}

TEST(FastgapDecoder, ThresholdCensoring) {
    std::mt19937_64 rng(5);
    auto dem = testutil::random_grid_dem(rng, 6, 5, 0.3);
    auto idx = testutil::make_index(dem);
    GapDecoder dec(idx, 1);
    for (int s = 0; s < 300; s++) {
        auto events = testutil::sample_dem(dem, rng);
        auto full = dec.decode(events, DecodeOptions());
        for (dist_int tau : {full.gap_lb / 2, full.gap_lb, full.gap_lb + 1, full.gap_lb * 2 + 1}) {
            DecodeOptions o;
            o.threshold = tau;
            auto r = dec.decode(events, o);
            if (r.censored) {
                EXPECT_EQ(r.gap_lb, tau);
                EXPECT_GE(full.gap_lb, tau);
            } else {
                EXPECT_EQ(r.gap_lb, full.gap_lb);
            }
            // Post-selection decision is never wrong.
            EXPECT_EQ(r.gap_lb >= tau, full.gap_lb >= tau);
        }
    }
}

TEST(FastgapDecoder, BatchAndThreadsMatchSingleShot) {
    auto dem = testutil::surface_code_dem(5, 5, 0.01, "rotated_memory_x");
    auto idx = testutil::make_index(dem);
    std::mt19937_64 rng(9);
    std::vector<std::vector<uint64_t>> shots;
    for (int s = 0; s < 300; s++)
        shots.push_back(testutil::sample_dem(dem, rng));
    // Many very large shots, so the default grain also takes the parallel search path.
    for (int s = 0; s < 40; s++) {
        std::vector<uint64_t> big;
        for (uint64_t d = s % 3; d < idx->num_detectors(); d += 3)
            big.push_back(d);
        if (big.size() % 2)
            big.pop_back();
        shots.push_back(big);
    }
    size_t bytes = (idx->num_detectors() + 7) / 8;
    auto packed = pack(shots, bytes);
    DecodeOptions o;
    o.upper_bound = true;
    o.check = false;
    o.keep_walk = true;
    GapDecoder single(idx, 1);
    std::vector<GapResult> batch1;
    single.decode_batch(packed.data(), shots.size(), bytes, o, batch1);
    // Every thread count and grain must reproduce the single-threaded search exactly, walk
    // included (ties are broken independently of the partition). Grain 1 forces the parallel
    // search even on shots with a handful of defects.
    for (size_t threads : {2u, 3u, 4u, 7u}) {
        for (size_t grain : {1u, 16u}) {
            GapDecoder multi(idx, threads);
            multi.parallel_grain = grain;
            std::vector<GapResult> batch;
            multi.decode_batch(packed.data(), shots.size(), bytes, o, batch);
            for (size_t s = 0; s < shots.size(); s++) {
                // Unsorted events exercise the target layout of each thread.
                auto shuffled = shots[s];
                std::shuffle(shuffled.begin(), shuffled.end(), rng);
                auto r = single.decode(shots[s], o);
                auto rm = multi.decode(shuffled, o);
                auto rs = single.decode(shuffled, o);
                for (auto* other : {&batch1[s], &batch[s]}) {
                    EXPECT_EQ(other->gap_lb, r.gap_lb);
                    EXPECT_EQ(other->gap_ub, r.gap_ub);
                    EXPECT_EQ(other->w_star, r.w_star);
                    EXPECT_EQ(other->prediction, r.prediction);
                    EXPECT_EQ(other->walk_simple, r.walk_simple);
                    ASSERT_EQ(other->walk.size(), r.walk.size());
                    for (size_t h = 0; h < r.walk.size(); h++) {
                        EXPECT_EQ(other->walk[h].from, r.walk[h].from);
                        EXPECT_EQ(other->walk[h].to, r.walk[h].to);
                        EXPECT_EQ(other->walk[h].price, r.walk[h].price);
                    }
                }
                EXPECT_EQ(rm.gap_lb, rs.gap_lb);
                EXPECT_EQ(rm.gap_ub, rs.gap_ub);
                EXPECT_EQ(rm.prediction, rs.prediction);
                ASSERT_EQ(rm.walk.size(), rs.walk.size());
                for (size_t h = 0; h < rs.walk.size(); h++) {
                    EXPECT_EQ(rm.walk[h].from, rs.walk[h].from);
                    EXPECT_EQ(rm.walk[h].to, rs.walk[h].to);
                }
                // Early termination agrees too.
                DecodeOptions ot = o;
                ot.threshold = r.gap_lb / 2 + 1;
                auto c1 = single.decode(shots[s], ot);
                auto cm = multi.decode(shots[s], ot);
                EXPECT_EQ(c1.gap_lb, cm.gap_lb);
                EXPECT_EQ(c1.censored, cm.censored);
            }
        }
    }
}

TEST(FastgapDecoder, SurfaceCodeSoundness) {
    for (auto task : {"rotated_memory_x", "rotated_memory_z"}) {
        for (uint32_t d : {3u, 5u}) {
            auto dem = testutil::surface_code_dem(d, d, 0.01, task);
            auto idx = testutil::make_index(dem);
            GapDecoder dec(idx, 1);
            ExactBaseline baseline(idx);
            std::mt19937_64 rng(d);
            DecodeOptions o;
            o.upper_bound = true;
            o.check = true;
            size_t blossoms = 0;
            for (int s = 0; s < 300; s++) {
                auto events = testutil::sample_dem(dem, rng);
                auto r = dec.decode(events, o);
                auto e = baseline.decode(events);
                ASSERT_LE(r.gap_lb, e.gap);
                if (r.walk_simple)
                    ASSERT_GE(r.gap_ub, e.gap);
                ASSERT_EQ(r.w_star, e.w_star);
                blossoms += r.num_blossoms > 0;
            }
            EXPECT_GT(blossoms, 0u);
        }
    }
}

