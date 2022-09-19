#include "pymatching/fill_match/driver/stim_io.h"

#include "gtest/gtest.h"

TEST(StimIO, ProbabilityGraph) {
    pm::ProbabilityGraph g(10, 64);
    g.add_or_merge_edge(0, 1, 0.2, {0});
    g.add_or_merge_edge(1, 2, 0.1, {0, 1});
    g.add_or_merge_edge(1, 2, 0.05, {0, 1});
    ASSERT_EQ(g.nodes[0][0].node, &g.nodes[1]);
    ASSERT_EQ(g.nodes[1][0].node, &g.nodes[0]);
    ASSERT_EQ(g.nodes[1][1].node, &g.nodes[2]);
    ASSERT_EQ(g.nodes[2][0].node, &g.nodes[1]);
    ASSERT_EQ(g.nodes[1][1].probability, 0.14);
    ASSERT_EQ(g.nodes[2][0].probability, 0.14);
    g.add_or_merge_boundary_edge(0, 0.3, {0, 2});
    ASSERT_EQ(g.nodes[0][1].probability, 0.3);
    std::vector<size_t> obs_expected = {0, 2};
    ASSERT_EQ(g.nodes[0][1].observables, obs_expected);
    g.add_or_merge_boundary_edge(0, 0.4, {0});
    ASSERT_FLOAT_EQ(g.nodes[0][1].probability, 0.46);
    ASSERT_EQ(g.nodes[0][1].observables, obs_expected);
}
