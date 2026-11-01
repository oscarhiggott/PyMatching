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

#include "pymatching/fastgap/gap_index.h"

#include <gtest/gtest.h>

#include <cstdio>

#include "pymatching/fastgap/fastgap_test_util.test.h"
#include "pymatching/fastgap/gap_decoder.h"

using namespace pm::fastgap;

TEST(FastgapIndex, DHatIsALowerBound) {
    std::mt19937_64 rng(42);
    auto dem = testutil::random_grid_dem(rng, 6, 5, 0.3);
    auto idx = testutil::make_index(dem);
    for (uint32_t u = 0; u < idx->graph.num_nodes; u++) {
        auto row = idx->graph.dijkstra({{u, 0}});
        for (uint32_t v = 0; v < idx->graph.num_nodes; v++) {
            bool exact;
            dist_int dh = idx->d_hat(u, v, exact);
            ASSERT_LE(dh, row[v]);
            ASSERT_EQ(dh, idx->lower_bound(u, v));
            ASSERT_EQ(exact, u == v);
        }
    }
}

TEST(FastgapIndex, LandmarkSetsAreLowerBounds) {
    auto dem = testutil::surface_code_dem(5, 5, 0.005, "rotated_memory_z");
    auto idx = GapIndex::create_from_dem(dem.str(), 0);
    std::vector<LandmarkSource> sources = {{"first_half", {}}, {"single", {7}}};
    for (uint32_t v = 0; v < idx->graph.num_nodes / 2; v++)
        sources[0].nodes.push_back(v);
    idx->set_landmarks(sources, "{}");
    for (uint32_t u = 0; u < idx->graph.num_nodes; u += 3) {
        auto row = idx->graph.dijkstra({{u, 0}});
        for (uint32_t v = 0; v < idx->graph.num_nodes; v++)
            ASSERT_LE(idx->lower_bound(u, v), row[v]);
    }
}

TEST(FastgapIndex, SaveLoadRoundTripIsBitIdentical) {
    auto dem = testutil::surface_code_dem(3, 3, 0.01, "rotated_memory_x");
    auto idx = testutil::make_index(dem);
    std::string path = testing::TempDir() + "fastgap_roundtrip.idx";
    idx->save(path);
    auto loaded = GapIndex::load(path);
    EXPECT_EQ(loaded->potentials, idx->potentials);
    EXPECT_EQ(loaded->strategy_hash, idx->strategy_hash);
    EXPECT_EQ(loaded->strategy_config, idx->strategy_config);
    EXPECT_EQ(loaded->d_lr, idx->d_lr);
    GapDecoder a(idx, 1), b(loaded, 1);
    std::mt19937_64 rng(3);
    DecodeOptions o;
    o.upper_bound = true;
    for (int s = 0; s < 200; s++) {
        auto events = testutil::sample_dem(dem, rng);
        auto ra = a.decode(events, o);
        auto rb = b.decode(events, o);
        ASSERT_EQ(ra.gap_lb, rb.gap_lb);
        ASSERT_EQ(ra.gap_ub, rb.gap_ub);
        ASSERT_EQ(ra.prediction, rb.prediction);
        ASSERT_EQ(ra.w_star, rb.w_star);
    }
    // A corrupted file is rejected.
    {
        FILE* f = fopen(path.c_str(), "r+b");
        ASSERT_NE(f, nullptr);
        fseek(f, 0, SEEK_SET);
        fputc('X', f);
        fclose(f);
    }
    EXPECT_THROW(GapIndex::load(path), std::runtime_error);
    std::remove(path.c_str());
}

TEST(FastgapIndex, ReferenceEdgesPath) {
    // Repetition code in the L/R model given as reference edges, with a detector-free L-R edge.
    std::vector<RefEdge> edges = {
        {0, SIZE_MAX, true, 2.0, -1},
        {0, 1, false, 2.0, -1},
        {1, SIZE_MAX, false, 2.0, -1},
        {SIZE_MAX, SIZE_MAX, false, 3.0, -1},
    };
    auto idx = GapIndex::create_from_reference_edges(3, 2, edges);
    idx->set_landmarks({}, "{}");
    // Bare chain through the interior: 6; direct edge: 3.
    EXPECT_EQ(idx->d_lr, (dist_int)std::llround(3.0 * idx->graph.normalising_constant));
    GapDecoder dec(idx, 1);
    auto r = dec.decode({2}, DecodeOptions());  // the obs detector is ignored
    EXPECT_EQ(r.num_defects, 0u);
    EXPECT_EQ(r.gap_lb, idx->d_lr);
}
