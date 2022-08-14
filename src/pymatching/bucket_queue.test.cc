#include "pymatching/bucket_queue.h"

#include <gtest/gtest.h>
#include <random>

TEST(bit_bucket_queue, basic_usage) {
    pm::bit_bucket_queue q;
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(q.cur_time, 0);

    pm::TentativeEvent out{};
    ASSERT_FALSE(q.try_pop_valid(&out));

    q.enqueue(pm::TentativeEvent(3));
    q.enqueue(pm::TentativeEvent(6, 0xDEAD));
    q.enqueue(pm::TentativeEvent(9));
    ASSERT_EQ(q.size(), 3);
    ASSERT_EQ(q.cur_time, 0);

    ASSERT_EQ(q.force_pop_valid(), pm::TentativeEvent(3, 0));
    ASSERT_EQ(q.size(), 2);
    ASSERT_EQ(q.cur_time, 3);

    ASSERT_EQ(q.force_pop_valid(), pm::TentativeEvent(9, 0));
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(q.cur_time, 9);
}

TEST(bit_bucket_queue, sorts_fuzz) {
    pm::bit_bucket_queue q;
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
        auto t = q.force_pop_valid().time;
        ASSERT_EQ(t, s[k]) << k;
    }

    pm::TentativeEvent out{};
    ASSERT_FALSE(q.try_pop_valid(&out));
}

TEST(bit_bucket_queue, bucket_for) {
    pm::bit_bucket_queue q;
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
