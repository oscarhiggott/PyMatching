#include "pymatching/bucket_queue.h"

#include <gtest/gtest.h>
#include <random>

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

TEST(bucket_queue, std_usage) {
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

TEST(bucket_queue, new_usage) {
    pm::bucket_queue<1> q(16);
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

TEST(bucket_queue, new2_usage) {
    pm::bucket_queue<2> q(16);
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

TEST(bucket_queue, new_sorts) {
    pm::bucket_queue<1> q(16);
    std::mt19937 rng(0); // NOLINT(cert-msc51-cpp)

    std::vector<pm::time_int> s;
    for (size_t k = 0; k < 1000; k++) {
        auto v = (pm::time_int)(rng() & 0x7FFFFFFF);
        s.push_back(v);
        q.enqueue(pm::TentativeEvent(v));
    }
    ASSERT_TRUE(q.satisfies_heap_invariant());
    std::sort(s.begin(), s.end());

    for (size_t k = 0; k < s.size(); k++) {
        auto t = q.force_pop().time;
        ASSERT_EQ(t, s[k]) << k;
    }

    pm::TentativeEvent out{};
    ASSERT_FALSE(q.try_pop(&out));
}


TEST(bucket_queue, new2_sorts) {
    pm::bucket_queue<2> q(16);
    std::mt19937 rng(0); // NOLINT(cert-msc51-cpp)

    std::vector<pm::time_int> s;
    for (size_t k = 0; k < 1000; k++) {
        auto v = (pm::time_int)(rng() & 0x7FFFFFFF);
        s.push_back(v);
        q.enqueue(pm::TentativeEvent(v));
    }
    std::sort(s.begin(), s.end());
    ASSERT_TRUE(q.satisfies_invariants());

    for (size_t k = 0; k < s.size(); k++) {
        auto t = q.force_pop().time;
        ASSERT_EQ(t, s[k]) << k;
    }

    pm::TentativeEvent out{};
    ASSERT_FALSE(q.try_pop(&out));
}
