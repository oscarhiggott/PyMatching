#include "pymatching/sparse_blossom/flooder/match.h"

#include <gtest/gtest.h>

using namespace pm;

TEST(MatchTest, equality) {
    Match m1{nullptr, {nullptr, nullptr, 5}};
    Match m2{nullptr, {nullptr, nullptr, 6}};
    ASSERT_TRUE(m1 != m2);
    ASSERT_TRUE(m1 == m1);
    ASSERT_FALSE(m1 == m2);
    ASSERT_FALSE(m1 != m1);
}
