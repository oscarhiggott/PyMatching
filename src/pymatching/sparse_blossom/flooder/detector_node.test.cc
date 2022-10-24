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

#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/flooder/graph.h"

using namespace pm;

TEST(DetectorNode, IndexOfNeighber) {
    DetectorNode d1;
    DetectorNode d2;
    d1.neighbors.push_back(nullptr);
    d1.neighbors.push_back(&d2);
    ASSERT_EQ(d1.index_of_neighbor(nullptr), 0);
    ASSERT_EQ(d1.index_of_neighbor(&d2), 1);
}

TEST(DetectorNode, compute_wrapped_radius_within_layer_at_time) {
    GraphFillRegion owner;
    GraphFillRegion above;
    GraphFillRegion top;
    DetectorNode node;
    node.reached_from_source = &node;
    node.region_that_arrived = &owner;
    node.region_that_arrived_top = &top;
    node.radius_of_arrival = 5;
    owner.blossom_parent = &above;
    above.blossom_parent = &top;
    top.blossom_parent = nullptr;
    owner.radius = VaryingCT::frozen(10);
    above.radius = VaryingCT::frozen(20);
    top.radius = VaryingCT::growing_value_at_time(0, 30);

    ASSERT_EQ(node.compute_wrapped_radius(), 25);

    ASSERT_EQ(node.compute_local_radius_at_time_bounded_by_region(25, owner), 5);
    ASSERT_EQ(node.compute_local_radius_at_time_bounded_by_region(30, owner), 5);
    ASSERT_EQ(node.compute_local_radius_at_time_bounded_by_region(40, owner), 5);

    ASSERT_EQ(node.compute_local_radius_at_time_bounded_by_region(25, above), 25);
    ASSERT_EQ(node.compute_local_radius_at_time_bounded_by_region(30, above), 25);
    ASSERT_EQ(node.compute_local_radius_at_time_bounded_by_region(40, above), 25);

    ASSERT_EQ(node.compute_local_radius_at_time_bounded_by_region(25, top), 20);
    ASSERT_EQ(node.compute_local_radius_at_time_bounded_by_region(30, top), 25);
    ASSERT_EQ(node.compute_local_radius_at_time_bounded_by_region(40, top), 35);
}

TEST(DetectorNode, compute_stitch_radius_at_time_bounded_by_region_towards_neighbor_across_collision) {
    GraphFillRegion left;
    GraphFillRegion right;
    GraphFillRegion parent;
    DetectorNode left_node;
    DetectorNode right_node;
    left_node.reached_from_source = &left_node;
    right_node.reached_from_source = &right_node;
    left_node.region_that_arrived = &left;
    right_node.region_that_arrived = &right;
    left.blossom_parent = &parent;
    right.blossom_parent = &parent;
    parent.blossom_parent = nullptr;
    left.radius = VaryingCT::frozen(5);
    right.radius = VaryingCT::frozen(5);
    parent.radius = VaryingCT::growing_value_at_time(0, 5);
    left_node.neighbors.push_back(&right_node);
    right_node.neighbors.push_back(&left_node);
    left_node.neighbor_weights.push_back(20);
    right_node.neighbor_weights.push_back(20);
    left_node.radius_of_arrival = 0;
    right_node.radius_of_arrival = 0;

    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(5, left, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(6, left, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(9, left, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(10, left, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(11, left, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(100, left, 0), 5);

    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(5, parent, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(6, parent, 0), 6);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(9, parent, 0), 9);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(10, parent, 0), 10);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(11, parent, 0), 10);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(102, parent, 0), 10);
    ASSERT_EQ(right_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(102, parent, 0), 10);
}

TEST(DetectorNode, compute_stitch_radius_at_time_bounded_by_region_towards_neighbor_across_skewed_collision) {
    GraphFillRegion left;
    GraphFillRegion right;
    GraphFillRegion parent;
    DetectorNode left_node;
    DetectorNode right_node;
    left_node.reached_from_source = &left_node;
    right_node.reached_from_source = &right_node;
    left_node.region_that_arrived = &left;
    right_node.region_that_arrived = &right;
    left.blossom_parent = &parent;
    right.blossom_parent = &parent;
    parent.blossom_parent = nullptr;
    left.radius = VaryingCT::frozen(5);
    right.radius = VaryingCT::frozen(8);
    parent.radius = VaryingCT::growing_value_at_time(0, 5);
    left_node.neighbors.push_back(&right_node);
    right_node.neighbors.push_back(&left_node);
    left_node.neighbor_weights.push_back(20);
    right_node.neighbor_weights.push_back(20);
    left_node.radius_of_arrival = 0;
    right_node.radius_of_arrival = 0;

    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(5, left, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(6, left, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(9, left, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(10, left, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(11, left, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(100, left, 0), 5);

    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(5, parent, 0), 5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(6, parent, 0), 6);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(8, parent, 0), 8);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(9, parent, 0), 8.5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(10, parent, 0), 8.5);
    ASSERT_EQ(left_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(100, parent, 0), 8.5);
    ASSERT_EQ(right_node.compute_stitch_radius_at_time_bounded_by_region_towards_neighbor(100, parent, 0), 11.5);
}
