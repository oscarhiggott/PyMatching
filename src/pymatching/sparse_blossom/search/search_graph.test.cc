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

#include "pymatching/sparse_blossom/search/search_graph.h"

#include <gtest/gtest.h>
#include <vector>

TEST(SearchGraph, AddEdge) {
    pm::SearchGraph g(4);
    g.add_edge(0, 1, 2, {0});
    g.add_edge(1, 2, 3, {0, 2});
    g.add_edge(0, 3, 10, {1, 3});
    ASSERT_EQ(g.nodes[0].neighbors[0], &g.nodes[1]);
    ASSERT_EQ(g.nodes[1].neighbors[0], &g.nodes[0]);
    ASSERT_EQ(g.nodes[3].neighbors[0], &g.nodes[0]);
    ASSERT_EQ(g.nodes[0].neighbors[1], &g.nodes[3]);
    ASSERT_EQ(g.nodes[0].neighbor_weights[0], 2);
    ASSERT_EQ(g.nodes[0].neighbor_weights[1], 10);
    std::vector<size_t> v1 = {0};
    ASSERT_EQ(g.nodes[1].neighbor_observable_indices[0], v1);
    std::vector<size_t> v2 = {1, 3};
    ASSERT_EQ(g.nodes[3].neighbor_observable_indices[0], v2);
}

TEST(SearchGraph, AddBoundaryEdge) {
    pm::SearchGraph g(6);
    g.add_edge(0, 1, 2, {0, 1});
    g.add_boundary_edge(0, 7, {2});
    g.add_boundary_edge(5, 10, {0, 1, 3});
    ASSERT_EQ(g.nodes[0].neighbors[0], nullptr);
    ASSERT_EQ(g.nodes[0].neighbors[1], &g.nodes[1]);
    ASSERT_EQ(g.nodes[5].neighbors[0], nullptr);
    std::vector<size_t> v1 = {2};
    ASSERT_EQ(g.nodes[0].neighbor_observable_indices[0], v1);
    ASSERT_EQ(g.nodes[0].neighbor_weights[0], 7);
}
