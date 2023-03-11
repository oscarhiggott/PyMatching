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

#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"

#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/flooder/graph_fill_region.test.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/region_edge.h"

using namespace pm;

TEST(GraphFillRegion, BlossomBuilder) {
    GraphFillTestData x(4);

    auto t1 = x.b(-1, -1, {x.b(1, 2, {}, false), x.b(2, 3, {}, false), x.b(3, 1, {}, false)}, true);
    auto t2 = x.b(-1, -1, {x.b(1, 2, {}, false), x.b(2, 3, {}, false), x.b(3, 1, {}, false)}, true);
    ASSERT_EQ(*t1.region, *t2.region);
}

TEST(GraphFillRegion, TopRegion) {
    GraphFillTestData x(10);
    auto t =
        x.b(-1,
            -1,
            {x.b(0, 1, {}, false),
             x.b(1, 2, {}, false),
             x.b(2, 0, {x.b(3, 4, {}, false), x.b(4, 5, {}, false), x.b(5, 3, {}, false)}, false)},
            true);

    ASSERT_EQ(t.region, t.region->blossom_children[2].region->blossom_children[0].region->blossom_parent_top);
    ASSERT_EQ(t.region, t.region->blossom_children[1].region->blossom_parent_top);
}

TEST(GraphFillRegion, AddMatch) {
    GraphFillRegion r1;
    GraphFillRegion r2;
    std::vector<DetectorNode> ds(2);
    r1.add_match(&r2, CompressedEdge{&ds[0], &ds[1], 5});
    ASSERT_EQ(r1.match.region, &r2);
    ASSERT_EQ(r1.match.edge, r2.match.edge.reversed());
    ASSERT_EQ(r1.match.edge, (CompressedEdge{&ds[0], &ds[1], 5}));
}

TEST(GraphFillRegion, descendance_ordering_le) {
    std::array<GraphFillRegion, 6> r;
    r[0].blossom_parent = &r[2];
    r[1].blossom_parent = &r[2];
    r[2].blossom_parent = &r[3];
    r[3].blossom_parent = nullptr;
    r[4].blossom_parent = &r[5];
    r[5].blossom_parent = nullptr;

    ASSERT_TRUE(r[0] <= r[0]);
    ASSERT_FALSE(r[0] <= r[1]);
    ASSERT_TRUE(r[0] <= r[2]);
    ASSERT_TRUE(r[0] <= r[3]);
    ASSERT_FALSE(r[0] <= r[4]);
    ASSERT_FALSE(r[0] <= r[5]);

    ASSERT_FALSE(r[1] <= r[0]);
    ASSERT_TRUE(r[1] <= r[1]);
    ASSERT_TRUE(r[1] <= r[2]);
    ASSERT_TRUE(r[1] <= r[3]);
    ASSERT_FALSE(r[1] <= r[4]);
    ASSERT_FALSE(r[1] <= r[5]);

    ASSERT_FALSE(r[2] <= r[0]);
    ASSERT_FALSE(r[2] <= r[1]);
    ASSERT_TRUE(r[2] <= r[2]);
    ASSERT_TRUE(r[2] <= r[3]);
    ASSERT_FALSE(r[2] <= r[4]);
    ASSERT_FALSE(r[2] <= r[5]);

    ASSERT_FALSE(r[3] <= r[0]);
    ASSERT_FALSE(r[3] <= r[1]);
    ASSERT_FALSE(r[3] <= r[2]);
    ASSERT_TRUE(r[3] <= r[3]);
    ASSERT_FALSE(r[3] <= r[4]);
    ASSERT_FALSE(r[3] <= r[5]);

    ASSERT_FALSE(r[4] <= r[0]);
    ASSERT_FALSE(r[4] <= r[1]);
    ASSERT_FALSE(r[4] <= r[2]);
    ASSERT_FALSE(r[4] <= r[3]);
    ASSERT_TRUE(r[4] <= r[4]);
    ASSERT_TRUE(r[4] <= r[5]);

    ASSERT_FALSE(r[5] <= r[0]);
    ASSERT_FALSE(r[5] <= r[1]);
    ASSERT_FALSE(r[5] <= r[2]);
    ASSERT_FALSE(r[5] <= r[3]);
    ASSERT_FALSE(r[5] <= r[4]);
    ASSERT_TRUE(r[5] <= r[5]);
}

TEST(GraphFillRegion, descendance_ordering_lt) {
    std::array<GraphFillRegion, 6> r;
    r[0].blossom_parent = &r[2];
    r[1].blossom_parent = &r[2];
    r[2].blossom_parent = &r[3];
    r[3].blossom_parent = nullptr;
    r[4].blossom_parent = &r[5];
    r[5].blossom_parent = nullptr;

    ASSERT_FALSE(r[0] < r[0]);
    ASSERT_FALSE(r[0] < r[1]);
    ASSERT_TRUE(r[0] < r[2]);
    ASSERT_TRUE(r[0] < r[3]);
    ASSERT_FALSE(r[0] < r[4]);
    ASSERT_FALSE(r[0] < r[5]);

    ASSERT_FALSE(r[1] < r[0]);
    ASSERT_FALSE(r[1] < r[1]);
    ASSERT_TRUE(r[1] < r[2]);
    ASSERT_TRUE(r[1] < r[3]);
    ASSERT_FALSE(r[1] < r[4]);
    ASSERT_FALSE(r[1] < r[5]);

    ASSERT_FALSE(r[2] < r[0]);
    ASSERT_FALSE(r[2] < r[1]);
    ASSERT_FALSE(r[2] < r[2]);
    ASSERT_TRUE(r[2] < r[3]);
    ASSERT_FALSE(r[2] < r[4]);
    ASSERT_FALSE(r[2] < r[5]);

    ASSERT_FALSE(r[3] < r[0]);
    ASSERT_FALSE(r[3] < r[1]);
    ASSERT_FALSE(r[3] < r[2]);
    ASSERT_FALSE(r[3] < r[3]);
    ASSERT_FALSE(r[3] < r[4]);
    ASSERT_FALSE(r[3] < r[5]);

    ASSERT_FALSE(r[4] < r[0]);
    ASSERT_FALSE(r[4] < r[1]);
    ASSERT_FALSE(r[4] < r[2]);
    ASSERT_FALSE(r[4] < r[3]);
    ASSERT_FALSE(r[4] < r[4]);
    ASSERT_TRUE(r[4] < r[5]);

    ASSERT_FALSE(r[5] < r[0]);
    ASSERT_FALSE(r[5] < r[1]);
    ASSERT_FALSE(r[5] < r[2]);
    ASSERT_FALSE(r[5] < r[3]);
    ASSERT_FALSE(r[5] < r[4]);
    ASSERT_FALSE(r[5] < r[5]);
}
