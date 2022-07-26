#include <gtest/gtest.h>
#include "compressed_edge.h"
#include "graph.h"

TEST(CompressedEdge, Reversed) {
    pm::Graph g(2);
    auto x = pm::CompressedEdge(&g.nodes[0], &g.nodes[1], 3);
    ASSERT_EQ(x.reversed(), pm::CompressedEdge(&g.nodes[1], &g.nodes[0], 3));
}

TEST(CompressedEdge, MergedWith) {
    pm::Graph g(3);
    auto x = pm::CompressedEdge(&g.nodes[0], &g.nodes[1], 5);
    auto y = pm::CompressedEdge(&g.nodes[1], &g.nodes[2], 6);
    ASSERT_EQ(x.merged_with(y), pm::CompressedEdge(&g.nodes[0], &g.nodes[2], 3));
}
