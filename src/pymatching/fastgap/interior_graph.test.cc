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

#include "pymatching/fastgap/interior_graph.h"

#include <gtest/gtest.h>

#include "pymatching/fastgap/fastgap_test_util.test.h"

using namespace pm::fastgap;

TEST(FastgapInteriorGraph, ReadsDecoderGraphWithSides) {
    stim::DetectorErrorModel dem(R"DEM(
        error(0.1) D0 L0
        error(0.1) D0 D1
        error(0.2) D1 D2
        error(0.1) D2
    )DEM");
    auto idx = GapIndex::create_from_dem(dem.str(), 0);
    auto& g = idx->graph;
    EXPECT_EQ(g.num_nodes, 3u);
    EXPECT_EQ(g.boundary_side[0], SIDE_L);
    EXPECT_EQ(g.boundary_side[1], SIDE_NONE);
    EXPECT_EQ(g.boundary_side[2], SIDE_R);
    // h_L(2) = w(0,b) + dist(0,2); d_LR is the full chain.
    auto row = g.dijkstra({{0, 0}});
    EXPECT_EQ(idx->h_l[2], g.boundary_weight[0] + row[2]);
    EXPECT_EQ(idx->h_l[0], g.boundary_weight[0]);
    EXPECT_EQ(idx->d_lr, idx->h_l[2] + g.boundary_weight[2]);
    EXPECT_EQ(idx->d_lr, idx->h_r[0] + g.boundary_weight[0]);
}

TEST(FastgapInteriorGraph, AStarMatchesDijkstra) {
    std::mt19937_64 rng(3);
    for (int trial = 0; trial < 5; trial++) {
        auto dem = testutil::random_grid_dem(rng, 7, 6);
        auto idx = testutil::make_index(dem, 3);
        InteriorGraph::AStarScratch scratch;
        for (uint32_t a = 0; a < idx->graph.num_nodes; a += 5) {
            auto row = idx->graph.dijkstra({{a, 0}});
            for (uint32_t b = 0; b < idx->graph.num_nodes; b++) {
                dist_int d = idx->graph.astar(b, a, [&](uint32_t x) { return idx->lower_bound(x, a); }, scratch);
                ASSERT_EQ(d, row[b]);
                ASSERT_LE(idx->lower_bound(a, b), row[b]);
            }
        }
    }
}

TEST(FastgapInteriorGraph, DistancesNeverTransitTheBoundary) {
    // Two detectors that are close via the boundary but far through the interior.
    stim::DetectorErrorModel dem(R"DEM(
        error(0.3) D0 L0
        error(0.001) D0 D1
        error(0.3) D1 L0
    )DEM");
    auto idx = GapIndex::create_from_dem(dem.str(), 0);
    auto row = idx->graph.dijkstra({{0, 0}});
    EXPECT_EQ(row[1], idx->graph.weights[idx->graph.offsets[0]]);
    EXPECT_GT(row[1], idx->graph.boundary_weight[0] + idx->graph.boundary_weight[1]);
}
