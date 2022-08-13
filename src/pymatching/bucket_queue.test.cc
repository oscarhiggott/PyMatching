#include "pymatching/bucket_queue.h"

#include <gtest/gtest.h>

TEST(bucket_queue, basic_usage) {
    pm::bucket_queue<4> q(16);
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(q.cur_time, 0);

    pm::TentativeEvent out{};
    ASSERT_FALSE(q.try_pop(&out));

    q.enqueue(pm::TentativeEvent(3));
    q.enqueue(pm::TentativeEvent(6, 0xDEAD));
    q.enqueue(pm::TentativeEvent(9));
    ASSERT_EQ(q.size(), 3);
    ASSERT_EQ(q.cur_time, 0);

    ASSERT_EQ(q.force_pop(), pm::TentativeEvent(3, 0));
    ASSERT_EQ(q.size(), 2);
    ASSERT_EQ(q.cur_time, 3);

    ASSERT_EQ(q.force_pop(), pm::TentativeEvent(9, 0));
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(q.cur_time, 9);
}

TEST(bucket_queue, fallback_usage) {
    pm::bucket_queue<0> q(16);
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(q.cur_time, 0);

    pm::TentativeEvent out{};
    ASSERT_FALSE(q.try_pop(&out));

    q.enqueue(pm::TentativeEvent(3));
    q.enqueue(pm::TentativeEvent(6, 0xDEAD));
    q.enqueue(pm::TentativeEvent(9));
    ASSERT_EQ(q.size(), 3);
    ASSERT_EQ(q.cur_time, 0);

    ASSERT_EQ(q.force_pop(), pm::TentativeEvent(3, 0));
    ASSERT_EQ(q.size(), 2);
    ASSERT_EQ(q.cur_time, 3);

    ASSERT_EQ(q.force_pop(), pm::TentativeEvent(9, 0));
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(q.cur_time, 9);
}
