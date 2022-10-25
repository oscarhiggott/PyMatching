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

#include "pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.h"

#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/flooder/graph.h"

using namespace pm;

TEST(CompressedEdge, Reversed) {
    MatchingGraph g(2, 64);
    auto x = CompressedEdge{&g.nodes[0], &g.nodes[1], 3};
    ASSERT_EQ(x.reversed(), (CompressedEdge{&g.nodes[1], &g.nodes[0], 3}));
}

TEST(CompressedEdge, MergedWith) {
    MatchingGraph g(3, 64);
    auto x = CompressedEdge{&g.nodes[0], &g.nodes[1], 5};
    auto y = CompressedEdge{&g.nodes[1], &g.nodes[2], 6};
    ASSERT_EQ(x.merged_with(y), (CompressedEdge{&g.nodes[0], &g.nodes[2], 3}));
}
