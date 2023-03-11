// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"

#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/ints.h"

using namespace pm;

TEST(Varying, GrowingZeroRadiusAtTime) {
    auto x = VaryingCT::growing_varying_with_zero_distance_at_time(3);
    ASSERT_EQ(x.get_distance_at_time(3), 0);
    ASSERT_EQ(x.get_distance_at_time(10), 7);
    ASSERT_EQ(x.get_distance_at_time(0), -3);
}

TEST(Varying, GetDistanceAtTime) {
    auto x1 = VaryingCT((-10 << 2) + 1);
    ASSERT_EQ(x1.get_distance_at_time(15), 5);
    ASSERT_EQ(x1.get_distance_at_time(5), -5);
    ASSERT_EQ(x1.get_growing_distance_at_time(5), -5);
    auto x2 = VaryingCT((20 << 2) + 2);
    ASSERT_EQ(x2.get_distance_at_time(30), -10);
    ASSERT_EQ(x2.get_distance_at_time(10), 10);
    ASSERT_EQ(x2.get_shrinking_distance_at_time(10), 10);
    auto x3 = VaryingCT((120 << 2));
    ASSERT_EQ(x3.get_distance_at_time(1000), 120);
    ASSERT_EQ(x3.get_distance_at_time(30), 120);
    ASSERT_EQ(x3.y_intercept(), 120);
}

TEST(Varying, VaryingString) {
    auto x = Varying<int32_t>((12 << 2) + 2);
    ASSERT_EQ(x.str(), "12 - t");
}

TEST(Varying, ThenSlopeAt) {
    auto x = Varying<int32_t>::growing_varying_with_zero_distance_at_time(20);
    auto y = x.then_shrinking_at_time(30);
    ASSERT_EQ(y, Varying<int32_t>((40 << 2) + 2));
    auto y2 = y.then_frozen_at_time(50);
    ASSERT_EQ(y2, Varying<int32_t>((-10 << 2)));
    auto y3 = y2.then_growing_at_time(60);
    ASSERT_EQ(y3.str(), "-70 + t");
}

TEST(Varying, AdditionAndSubtraction) {
    auto x = VaryingCT((100 << 2) + 2);
    ASSERT_EQ(x.str(), "100 - t");
    ASSERT_EQ((x + 9).str(), "109 - t");
    ASSERT_EQ((x - 19).str(), "81 - t");
    auto x2 = VaryingCT((1000 << 2) + 1);
    ASSERT_EQ((x2 + 100).str(), "1100 + t");
    ASSERT_EQ((x2 - 200).str(), "800 + t");
    auto x3 = VaryingCT((50 << 2));
    ASSERT_EQ((x3 + 1000).str(), "1050");
    ASSERT_EQ((x3 + 1000).str(), "1050");
}

TEST(Varying, InplaceAdditionAndSubtraction) {
    auto x = VaryingCT();
    ASSERT_TRUE(x.is_frozen());
    ASSERT_EQ(x.y_intercept(), 0);

    x += 5;
    ASSERT_EQ(x.str(), "5");

    x -= 3;
    ASSERT_EQ(x.str(), "2");

    x = x.then_growing_at_time(0);
    ASSERT_EQ(x.str(), "2 + t");

    x += 5;
    ASSERT_EQ(x.str(), "7 + t");

    x -= 11;
    ASSERT_EQ(x.str(), "-4 + t");
}

TEST(Varying, from_base_and_growth) {
    ASSERT_EQ(VaryingCT::from_base_and_growth(5, -1).str(), "5 - t");
    ASSERT_EQ(VaryingCT::from_base_and_growth(5, 0).str(), "5");
    ASSERT_EQ(VaryingCT::from_base_and_growth(5, +1).str(), "5 + t");
}

TEST(Varying, growing_value_at_time) {
    auto x = VaryingCT::growing_value_at_time(10, 5);
    ASSERT_EQ(x.get_distance_at_time(5), 10);
    ASSERT_EQ(x.get_distance_at_time(6), 11);
    ASSERT_TRUE(x.is_growing());
}

TEST(Varying, shrinking_value_at_time) {
    auto x = Varying32::shrinking_value_at_time(11, 6);
    ASSERT_EQ(x.get_distance_at_time(6), 11);
    ASSERT_EQ(x.get_distance_at_time(7), 10);
    ASSERT_TRUE(x.is_shrinking());
}

TEST(Varying, frozen) {
    auto x = Varying32::frozen(11);
    ASSERT_EQ(x.get_distance_at_time(5), 11);
    ASSERT_EQ(x.get_distance_at_time(6), 11);
    ASSERT_TRUE(x.is_frozen());
}

TEST(Varying, TimeToXInterceptWhenAddedTo) {
    auto x1 = Varying32((-10 << 2) + 1);
    auto x2 = Varying32((-20 << 2) + 1);
    int32_t t = x1.time_of_x_intercept_when_added_to(x2 - 100);
    ASSERT_EQ(t, 65);

    auto x3 = Varying32((10 << 2));
    int32_t t2 = x1.time_of_x_intercept_when_added_to(x3 - 100);
    ASSERT_EQ(t2, 100);
}

TEST(Varying, TimeToXIntercept) {
    ASSERT_EQ(Varying64((-90 << 2) + 1).time_of_x_intercept(), 90);
    ASSERT_EQ(Varying64((120 << 2) + 2).time_of_x_intercept(), 120);
}

TEST(Varying, Equality) {
    ASSERT_TRUE(Varying32((-12345 << 2) + 1) == Varying32((-12345 << 2) + 1));
    ASSERT_FALSE(Varying32((-12 << 2) + 1) == Varying32((-13 << 2) + 1));
    ASSERT_TRUE(Varying32((-12345 << 2) + 1) != Varying32((-12345 << 2) + 2));
    ASSERT_TRUE(Varying32((-12345 << 2) + 1) != Varying32((-12346 << 2) + 1));
}
