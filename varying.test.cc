#include <gtest/gtest.h>
#include "varying.h"

TEST(Varying, GrowingZeroRadiusAtTime) {
    pm::Varying32 x = pm::Varying32::growing_varying_with_zero_distance_at_time(3);
    ASSERT_EQ(x.get_distance_at_time(3), 0);
    ASSERT_EQ(x.get_distance_at_time(10), 7);
    ASSERT_EQ(x.get_distance_at_time(0), -3);
}

TEST(Varying, GetDistanceAtTime) {
    auto x1 = pm::Varying32((-10 << 2) + 1);
    ASSERT_EQ(x1.get_distance_at_time(15), 5);
    ASSERT_EQ(x1.get_distance_at_time(5), -5);
    ASSERT_EQ(x1.get_growing_distance_at_time(5), -5);
    auto x2 = pm::Varying32((20 << 2) + 2);
    ASSERT_EQ(x2.get_distance_at_time(30), -10);
    ASSERT_EQ(x2.get_distance_at_time(10), 10);
    ASSERT_EQ(x2.get_shrinking_distance_at_time(10), 10);
    auto x3 = pm::Varying32((120 << 2));
    ASSERT_EQ(x3.get_distance_at_time(1000), 120);
    ASSERT_EQ(x3.get_distance_at_time(30), 120);
    ASSERT_EQ(x3.y_intercept(), 120);
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
    auto x = pm::Varying32((100 << 2) + 2);
    ASSERT_EQ(x.str(), "100 - t");
    ASSERT_EQ((x + 9).str(), "109 - t");
    ASSERT_EQ((x - 19).str(), "81 - t");
    auto x2 = pm::Varying32((1000 << 2) + 1);
    ASSERT_EQ((x2 + 100).str(), "1100 + t");
    ASSERT_EQ((x2 - 200).str(), "800 + t");
    auto x3 = pm::Varying32((50 << 2));
    ASSERT_EQ((x3 + 1000).str(), "1050");
    ASSERT_EQ((x3 + 1000).str(), "1050");
}

TEST(Varying, TimeToXInterceptWhenAddedTo) {
    auto x1 = pm::Varying32((-10 << 2) + 1);
    auto x2 = pm::Varying32 ((-20 << 2) + 1);
    int32_t t = x1.time_of_x_intercept_when_added_to(x2 - 100);
    ASSERT_EQ(t, 65);

    auto x3 = pm::Varying32 ((10 << 2));
    int32_t t2 = x1.time_of_x_intercept_when_added_to(x3 - 100);
    ASSERT_EQ(t2, 100);
}

TEST(Varying, TimeToXIntercept) {
    ASSERT_EQ(pm::Varying64((-90 << 2) + 1).time_of_x_intercept(), 90);
    ASSERT_EQ(pm::Varying64((120 << 2) + 2).time_of_x_intercept(), 120);
}

TEST(Varying, Equality) {
    ASSERT_TRUE(pm::Varying32((-12345 << 2) + 1) == pm::Varying32((-12345 << 2) + 1));
    ASSERT_FALSE(pm::Varying32((-12 << 2) + 1) == pm::Varying32((-13 << 2) + 1));
    ASSERT_TRUE(pm::Varying32((-12345 << 2) + 1) != pm::Varying32((-12345 << 2) + 2));
    ASSERT_TRUE(pm::Varying32((-12345 << 2) + 1) != pm::Varying32((-12346 << 2) + 1));
}
