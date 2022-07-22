#include <gtest/gtest.h>
#include "varying.h"

TEST(Varying, GrowingZeroRadiusAtTime) {
    pm::Varying x = pm::Varying<pm::time_int>::growing_varying_with_zero_distance_at_time(3);
    ASSERT_EQ(x.get_distance_at_time(3), 0);
    ASSERT_EQ(x.get_distance_at_time(10), 7);
    ASSERT_EQ(x.get_distance_at_time(0), -3);
}

TEST(Varying, VaryingString) {
    auto x = pm::Varying<pm::time_int>((12 << 2) + 2);
    ASSERT_EQ(x.str(), "12 - t");
}

TEST(Varying, ThenSlopeAt){
    auto x = pm::Varying<pm::time_int>::growing_varying_with_zero_distance_at_time(20);
    auto y = x.then_shrinking_at_time(30);
    ASSERT_EQ(y, pm::Varying<pm::time_int>((40 << 2) + 2));
    auto y2 = y.then_frozen_at_time(50);
    ASSERT_EQ(y2, pm::Varying<pm::time_int>((-10 << 2)));
    auto y3 = y2.then_growing_at_time(60);
    ASSERT_EQ(y3.str(), "-70 + t");
}

TEST(Varying, AdditionAndSubtraction) {
    auto x = pm::Varying<pm::time_int>((100 << 2) + 2);
    ASSERT_EQ(x.str(), "100 - t");
    ASSERT_EQ((x + 9).str(), "109 - t");
    ASSERT_EQ((x - 19).str(), "81 - t");
    auto x2 = pm::Varying<pm::time_int>((1000 << 2) + 1);
    ASSERT_EQ((x2 + 100).str(), "1100 + t");
    ASSERT_EQ((x2 - 200).str(), "800 + t");
    auto x3 = pm::Varying<pm::time_int>((50 << 2));
    ASSERT_EQ((x3 + 1000).str(), "1050");
    ASSERT_EQ((x3 + 1000).str(), "1050");
}
