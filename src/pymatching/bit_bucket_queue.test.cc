#include "pymatching/bit_bucket_queue.h"

#include <gtest/gtest.h>
#include <random>

using namespace pm;

TEST(bit_bucket_queue, basic_usage) {
    bit_bucket_queue<true> q;
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(q.cur_time, 0);

    ASSERT_EQ(q.dequeue(), TentativeEvent(cyclic_time_int{0}));

    q.enqueue(TentativeEvent(cyclic_time_int{9}));
    q.enqueue(TentativeEvent(cyclic_time_int{3}));
    ASSERT_EQ(q.size(), 2);
    ASSERT_EQ(q.cur_time, 0);

    ASSERT_EQ(q.dequeue(), TentativeEvent(cyclic_time_int{3}));
    ASSERT_EQ(q.dequeue(), TentativeEvent(cyclic_time_int{9}));
    ASSERT_EQ(q.size(), 0);
    ASSERT_EQ(q.cur_time, 9);

    ASSERT_EQ(q.dequeue(), TentativeEvent(cyclic_time_int{0}));
}

TEST(bit_bucket_queue, sorts_fuzz) {
    bit_bucket_queue<true> q;
    std::mt19937 rng(0); // NOLINT(cert-msc51-cpp)

    std::vector<cyclic_time_int> s;
    for (size_t k = 0; k < 1000; k++) {
        auto v = cyclic_time_int{rng() & 0x00007FFF};
        s.push_back(v);
        q.enqueue(TentativeEvent(v));
    }
    std::sort(s.begin(), s.end());
    ASSERT_TRUE(q.satisfies_invariants());

    for (size_t k = 0; k < s.size(); k++) {
        auto t = q.dequeue().time;
        ASSERT_EQ(t, s[k]) << k;
    }

    ASSERT_EQ(q.dequeue(), TentativeEvent(cyclic_time_int{0}));
}

TEST(bit_bucket_queue, bucket_for) {
    bit_bucket_queue<true> q;
    q.cur_time = 17;
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{17}), 0);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{18}), 2);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{19}), 2);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{20}), 3);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{21}), 3);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{22}), 3);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{23}), 3);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{24}), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{25}), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{26}), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{27}), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{28}), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{29}), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{30}), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{31}), 4);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{32}), 6);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{33}), 6);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{63}), 6);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{64}), 7);
    ASSERT_EQ(q.cur_bit_bucket_for(cyclic_time_int{65}), 7);
}

TEST(bit_bucket_queue, wraparound) {
    bit_bucket_queue<true> q;
    q.cur_time = INT32_MAX;
    q.enqueue(TentativeEvent(cyclic_time_int{INT32_MIN}));
    q.enqueue(TentativeEvent(cyclic_time_int{INT32_MIN + 1}));
    ASSERT_EQ(q.dequeue().time, INT32_MIN);
    ASSERT_EQ(q.dequeue().time, INT32_MIN + 1);
}

TEST(bit_bucket_queue, QueuedEventTracker) {
    QueuedEventTracker tracker;
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_NOT_QUEUED);

    ASSERT_FALSE(tracker.decide_to_dequeue(cyclic_time_int{0}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_NOT_QUEUED);

    tracker.invalidate();
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_NOT_QUEUED);

    ASSERT_TRUE(tracker.decide_to_queue(cyclic_time_int{5}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_QUEUED);
    ASSERT_EQ(tracker.queued_time, 5);

    ASSERT_FALSE(tracker.decide_to_queue(cyclic_time_int{5}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_QUEUED);
    ASSERT_EQ(tracker.queued_time, 5);

    ASSERT_FALSE(tracker.decide_to_queue(cyclic_time_int{6}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_QUEUED);
    ASSERT_EQ(tracker.queued_time, 5);

    ASSERT_TRUE(tracker.decide_to_queue(cyclic_time_int{4}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_QUEUED);
    ASSERT_EQ(tracker.queued_time, 4);

    tracker.invalidate();
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_QUEUED_BUT_IGNORE);
    ASSERT_EQ(tracker.queued_time, 4);

    ASSERT_FALSE(tracker.decide_to_queue(cyclic_time_int{4}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_QUEUED);
    ASSERT_EQ(tracker.queued_time, 4);

    tracker.invalidate();
    ASSERT_TRUE(tracker.decide_to_queue(cyclic_time_int{7}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_QUEUED);
    ASSERT_EQ(tracker.queued_time, 7);

    tracker.invalidate();
    ASSERT_FALSE(tracker.decide_to_dequeue(cyclic_time_int{8}));
    ASSERT_FALSE(tracker.decide_to_dequeue(cyclic_time_int{4}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_QUEUED_BUT_IGNORE);
    ASSERT_EQ(tracker.queued_time, 7);

    ASSERT_FALSE(tracker.decide_to_dequeue(cyclic_time_int{7}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_NOT_QUEUED);

    ASSERT_TRUE(tracker.decide_to_queue(cyclic_time_int{10}));
    ASSERT_FALSE(tracker.decide_to_dequeue(cyclic_time_int{11}));
    ASSERT_FALSE(tracker.decide_to_dequeue(cyclic_time_int{9}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_QUEUED);
    ASSERT_EQ(tracker.queued_time, 10);

    ASSERT_TRUE(tracker.decide_to_dequeue(cyclic_time_int{10}));
    ASSERT_EQ(tracker.queued_state, EVENT_SOURCE_NOT_QUEUED);
}
