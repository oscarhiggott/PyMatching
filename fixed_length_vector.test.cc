#include <gtest/gtest.h>
#include "fixed_length_vector.h"

TEST(FixedLengthVector, PushBack) {
    pm::FixedLengthVector<int, 10> a;
    ASSERT_EQ(a.size(), 0);
    a.push_back(61);
    ASSERT_EQ(a.size(), 1);
    ASSERT_EQ(a[0], 61);
    a.push_back(7);
    ASSERT_EQ(a[1], 7);
    ASSERT_EQ(a.size(), 2);
    a[0] = 100;
    ASSERT_EQ(a[0], 100);
    ASSERT_EQ(a.size(), 2);

    pm::Vec16<int> b;
    ASSERT_EQ(b.size(), 0);
    b.push_back(10);
    ASSERT_EQ(b.size(), 1);
    ASSERT_EQ(b[0], 10);
}
