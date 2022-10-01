#include "gtest/gtest.h"

#include "pymatching/fill_match/driver/io.h"

TEST(StimIO, IntermediateWeightedGraph) {
    pm::IntermediateWeightedGraph g(10, 64);
    g.add_or_merge_edge(0, 1, 0.2, {0});
    g.add_or_merge_edge(1, 2, 0.1, {0, 1});
    g.add_or_merge_edge(1, 2, 0.05, {0, 1});
    ASSERT_EQ(g.nodes[0][0].node, &g.nodes[1]);
    ASSERT_EQ(g.nodes[1][0].node, &g.nodes[0]);
    ASSERT_EQ(g.nodes[1][1].node, &g.nodes[2]);
    ASSERT_EQ(g.nodes[2][0].node, &g.nodes[1]);
    ASSERT_EQ(g.nodes[1][1].weight, 0.0024973997762458566);
    ASSERT_EQ(g.nodes[2][0].weight, 0.0024973997762458566);
    g.add_or_merge_boundary_edge(0, 0.3, {0, 2});
    ASSERT_EQ(g.nodes[0][1].weight, 0.3);
    std::vector<size_t> obs_expected = {0, 2};
    ASSERT_EQ(g.nodes[0][1].observables, obs_expected);
    g.add_or_merge_boundary_edge(0, 0.4, {0});
    ASSERT_FLOAT_EQ(g.nodes[0][1].weight, 0.05878938881188689);
    ASSERT_EQ(g.nodes[0][1].observables, obs_expected);
}

double merge_weights_via_probabilities(double a, double b){
    double p_a = 1/(1 + std::exp(a));
    double p_b = 1/(1 + std::exp(b));
    double p_both = p_a * (1 - p_b) + p_b * (1 - p_a);
    return std::log((1-p_both)/p_both);
}

TEST(StimIO, MergeWeights) {
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(10, 11), pm::merge_weights(10, 11));
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(-4.1, 5), pm::merge_weights(-4.1, 5));
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(25, -24), pm::merge_weights(25, -24));
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(-1, -2), pm::merge_weights(-1, -2));
    ASSERT_FLOAT_EQ(pm::merge_weights(1000, 0), 0);
}
