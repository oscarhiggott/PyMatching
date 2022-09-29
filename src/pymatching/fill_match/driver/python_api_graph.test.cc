#include <gtest/gtest.h>

#include <cmath>

#include "pymatching/fill_match/driver/python_api_graph.h"

double merge_weights_via_probabilities(double a, double b){
    double p_a = 1/(1 + std::exp(a));
    double p_b = 1/(1 + std::exp(b));
    double p_both = p_a * (1 - p_b) + p_b * (1 - p_a);
    return std::log((1-p_both)/p_both);
}

TEST(PythonApiGraph, MergeWeights) {
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(10, 11), pm::merge_weights(10, 11));
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(-4.1, 5), pm::merge_weights(-4.1, 5));
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(25, -24), pm::merge_weights(25, -24));
    ASSERT_FLOAT_EQ(merge_weights_via_probabilities(-1, -2), pm::merge_weights(-1, -2));
    ASSERT_FLOAT_EQ(pm::merge_weights(1000, 0), 0);
}
