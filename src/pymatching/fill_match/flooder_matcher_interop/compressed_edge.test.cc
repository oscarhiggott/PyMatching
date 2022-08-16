#include "pymatching/fill_match/flooder_matcher_interop/compressed_edge.h"

#include <gtest/gtest.h>

#include "pymatching/fill_match/flooder/graph.h"

using namespace pm;

TEST(CompressedEdge, Reversed) {
    MatchingGraph g(2);
    auto x = CompressedEdge{&g.nodes[0], &g.nodes[1], 3};
    ASSERT_EQ(x.reversed(), (CompressedEdge{&g.nodes[1], &g.nodes[0], 3}));
}

TEST(CompressedEdge, MergedWith) {
    MatchingGraph g(3);
    auto x = CompressedEdge{&g.nodes[0], &g.nodes[1], 5};
    auto y = CompressedEdge{&g.nodes[1], &g.nodes[2], 6};
    ASSERT_EQ(x.merged_with(y), (CompressedEdge{&g.nodes[0], &g.nodes[2], 3}));
}
