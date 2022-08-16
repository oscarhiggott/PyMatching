#include "pymatching/fill_match/flooder/graph.h"

#include <gtest/gtest.h>

#include "pymatching/fill_match/flooder/graph_fill_region.test.h"

using namespace pm;

TEST(DetectorNode, IndexOfNeighber) {
    DetectorNode d1;
    DetectorNode d2;
    d1.neighbors.push_back(nullptr);
    d1.neighbors.push_back(&d2);
    ASSERT_EQ(d1.index_of_neighbor(nullptr), 0);
    ASSERT_EQ(d1.index_of_neighbor(&d2), 1);
}
