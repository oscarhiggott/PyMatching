#include "pymatching/fill_match/flooder/graph.h"

#include <gtest/gtest.h>

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
    owner.radius = Varying32::from_base_and_growth(10, 0);
    above.radius = Varying32::from_base_and_growth(20, 0);
    top.radius = Varying32::from_point_and_slope(30, 0, +1);

    ASSERT_EQ(node.compute_wrapped_radius(), 25);

    ASSERT_EQ(node.compute_wrapped_radius_within_layer_at_time(&owner, 25), 5);
    ASSERT_EQ(node.compute_wrapped_radius_within_layer_at_time(&owner, 30), 5);
    ASSERT_EQ(node.compute_wrapped_radius_within_layer_at_time(&owner, 40), 5);

    ASSERT_EQ(node.compute_wrapped_radius_within_layer_at_time(&above, 25), 25);
    ASSERT_EQ(node.compute_wrapped_radius_within_layer_at_time(&above, 30), 25);
    ASSERT_EQ(node.compute_wrapped_radius_within_layer_at_time(&above, 40), 25);

    ASSERT_EQ(node.compute_wrapped_radius_within_layer_at_time(&top, 25), 20);
    ASSERT_EQ(node.compute_wrapped_radius_within_layer_at_time(&top, 30), 25);
    ASSERT_EQ(node.compute_wrapped_radius_within_layer_at_time(&top, 40), 35);
}
