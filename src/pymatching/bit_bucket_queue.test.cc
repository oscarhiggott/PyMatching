#include "pymatching/bit_bucket_queue.h"

#include <gtest/gtest.h>
#include <random>

TEST(bit_bucket_queue, basic_usage) {
    pm::bit_bucket_queue<true> q;
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(q.cur_time, 0);

    ASSERT_EQ(q.dequeue(), pm::TentativeEvent(0));

    q.enqueue(pm::TentativeEvent(9));
    q.enqueue(pm::TentativeEvent(3));
    ASSERT_EQ(q.size(), 2);
    ASSERT_EQ(q.cur_time, 0);

    ASSERT_EQ(q.dequeue(), pm::TentativeEvent(3));
    ASSERT_EQ(q.dequeue(), pm::TentativeEvent(9));
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(q.cur_time, 9);

    ASSERT_EQ(q.dequeue(), pm::TentativeEvent(0));
}

TEST(bit_bucket_queue, sorts_fuzz) {
    pm::bit_bucket_queue<true> q;
    std::mt19937 rng(0); // NOLINT(cert-msc51-cpp)

    std::vector<pm::time_int> s;
    for (size_t k = 0; k < 1000; k++) {
        auto v = (pm::time_int)(rng() & 0x3FFFFFFF);
        s.push_back(v);
        q.enqueue(pm::TentativeEvent(v));
    }
    std::sort(s.begin(), s.end());
    ASSERT_TRUE(q.satisfies_invariants());

    for (size_t k = 0; k < s.size(); k++) {
        auto t = q.dequeue().time;
        ASSERT_EQ(t, s[k]) << k;
    }

    ASSERT_EQ(q.dequeue(), pm::TentativeEvent(0));
}

TEST(bit_bucket_queue, bucket_for) {
    pm::bit_bucket_queue<true> q;
    q.cur_time = 17;
    ASSERT_EQ(q.cur_bit_bucket_for(17), 0);
    ASSERT_EQ(q.cur_bit_bucket_for(18), 2);
    ASSERT_EQ(q.cur_bit_bucket_for(19), 2);
    ASSERT_EQ(q.cur_bit_bucket_for(20), 3);
    ASSERT_EQ(q.cur_bit_bucket_for(21), 3);
    ASSERT_EQ(q.cur_bit_bucket_for(22), 3);
    ASSERT_EQ(q.cur_bit_bucket_for(23), 3);
    ASSERT_EQ(q.cur_bit_bucket_for(24), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(25), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(26), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(27), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(28), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(29), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(30), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(31), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(32), 6);
    ASSERT_EQ(q.cur_bit_bucket_for(33), 6);
    ASSERT_EQ(q.cur_bit_bucket_for(63), 6);
    ASSERT_EQ(q.cur_bit_bucket_for(64), 7);
    ASSERT_EQ(q.cur_bit_bucket_for(65), 7);
}

TEST(bit_bucket_queue, wraparound) {
    pm::bit_bucket_queue<true> q;
    q.cur_time = INT32_MAX;
    q.enqueue(pm::TentativeEvent(INT32_MIN));
    q.enqueue(pm::TentativeEvent(INT32_MIN + 1));
    ASSERT_EQ(q.dequeue().time, INT32_MIN);
    ASSERT_EQ(q.dequeue().time, INT32_MIN + 1);
}

TEST(bit_bucket_queue, cycle_comparison) {
    ASSERT_FALSE(pm::is_time_x_cyclebefore_y(1, 0));
    ASSERT_FALSE(pm::is_time_x_cycleatmost_y(1, 0));
    ASSERT_TRUE(pm::is_time_x_cycleatleast_y(1, 0));

    ASSERT_TRUE(pm::is_time_x_cyclebefore_y(0, 1));
    ASSERT_TRUE(pm::is_time_x_cycleatmost_y(0, 1));
    ASSERT_FALSE(pm::is_time_x_cycleatleast_y(0, 1));

    ASSERT_FALSE(pm::is_time_x_cyclebefore_y(0, 0));
    ASSERT_TRUE(pm::is_time_x_cycleatmost_y(0, 0));
    ASSERT_TRUE(pm::is_time_x_cycleatleast_y(0, 0));

    ASSERT_FALSE(pm::is_time_x_cyclebefore_y(INT32_MIN, INT32_MAX));
    ASSERT_FALSE(pm::is_time_x_cycleatmost_y(INT32_MIN, INT32_MAX));
    ASSERT_TRUE(pm::is_time_x_cycleatleast_y(INT32_MIN, INT32_MAX));

    ASSERT_TRUE(pm::is_time_x_cyclebefore_y(INT32_MAX, INT32_MIN));
    ASSERT_TRUE(pm::is_time_x_cycleatmost_y(INT32_MAX, INT32_MIN));
    ASSERT_FALSE(pm::is_time_x_cycleatleast_y(INT32_MAX, INT32_MIN));

    ASSERT_FALSE(pm::is_time_x_cyclebefore_y(INT32_MAX, INT32_MAX));
    ASSERT_TRUE(pm::is_time_x_cycleatmost_y(INT32_MAX, INT32_MAX));
    ASSERT_TRUE(pm::is_time_x_cycleatleast_y(INT32_MAX, INT32_MAX));
}
