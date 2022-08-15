#include "pymatching/fill_match/tracker/bit_bucket_queue.h"

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

TEST(bit_bucket_queue, QueuedEventTracker) {
    auto ev = [](int x){ return TentativeEvent{cyclic_time_int{x}}; };
    auto evs = [&](std::vector<int> x){
        std::vector<TentativeEvent> r;
        for (auto e : x) {
            r.push_back(ev(e));
        }
        return r;
    };
    bit_bucket_queue<true> queue;
    QueuedEventTracker tracker;

    ASSERT_EQ(tracker.has_queued_time, false);
    ASSERT_EQ(tracker.has_desired_time, false);

    tracker.set_desired_event(ev(5), queue);
    ASSERT_EQ(queue.cur_time, 0);
    ASSERT_EQ(queue.to_vector(), evs({5}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{5});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{5});

    tracker.set_desired_event(ev(6), queue);
    ASSERT_EQ(queue.cur_time, 0);
    ASSERT_EQ(queue.to_vector(), evs({5}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{5});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{6});

    tracker.set_desired_event(ev(4), queue);
    ASSERT_EQ(queue.cur_time, 0);
    ASSERT_EQ(queue.to_vector(), evs({4, 5}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{4});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{4});

    tracker.set_desired_event(ev(4), queue);
    ASSERT_EQ(queue.cur_time, 0);
    ASSERT_EQ(queue.to_vector(), evs({4, 5}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{4});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{4});

    tracker.set_desired_event(ev(3), queue);
    ASSERT_EQ(queue.cur_time, 0);
    ASSERT_EQ(queue.to_vector(), evs({3, 4, 5}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{3});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{3});

    auto deq = queue.dequeue();
    ASSERT_EQ(deq, ev(3));
    ASSERT_EQ(queue.cur_time, 3);
    ASSERT_EQ(queue.to_vector(), evs({4, 5}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{3});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{3});

    ASSERT_TRUE(tracker.dequeue_decision(deq, queue));
    ASSERT_EQ(queue.cur_time, 3);
    ASSERT_EQ(queue.to_vector(), evs({4, 5}));
    ASSERT_EQ(tracker.has_queued_time, false);
    ASSERT_EQ(tracker.has_desired_time, false);

    deq = queue.dequeue();
    ASSERT_EQ(deq, ev(4));
    ASSERT_EQ(queue.cur_time, 4);
    ASSERT_EQ(queue.to_vector(), evs({5}));
    ASSERT_EQ(tracker.has_queued_time, false);
    ASSERT_EQ(tracker.has_desired_time, false);

    ASSERT_FALSE(tracker.dequeue_decision(deq, queue));
    ASSERT_EQ(queue.cur_time, 4);
    ASSERT_EQ(queue.to_vector(), evs({5}));
    ASSERT_EQ(tracker.has_queued_time, false);
    ASSERT_EQ(tracker.has_desired_time, false);

    tracker.set_desired_event(ev(10), queue);
    ASSERT_EQ(queue.cur_time, 4);
    ASSERT_EQ(queue.to_vector(), evs({5, 10}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{10});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{10});

    deq = queue.dequeue();
    ASSERT_EQ(deq, ev(5));
    ASSERT_EQ(queue.cur_time, 5);
    ASSERT_EQ(queue.to_vector(), evs({10}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{10});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{10});

    ASSERT_FALSE(tracker.dequeue_decision(deq, queue));
    ASSERT_EQ(deq, ev(5));
    ASSERT_EQ(queue.cur_time, 5);
    ASSERT_EQ(queue.to_vector(), evs({10}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{10});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{10});

    tracker.set_desired_event(ev(11), queue);
    ASSERT_EQ(queue.cur_time, 5);
    ASSERT_EQ(queue.to_vector(), evs({10}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{10});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{11});

    deq = queue.dequeue();
    ASSERT_EQ(deq, ev(10));
    ASSERT_EQ(queue.cur_time, 10);
    ASSERT_EQ(queue.to_vector(), evs({}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{10});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{11});

    ASSERT_FALSE(tracker.dequeue_decision(deq, queue));
    ASSERT_EQ(queue.cur_time, 10);
    ASSERT_EQ(queue.to_vector(), evs({11}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{11});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{11});

    tracker.set_no_desired_event();
    ASSERT_EQ(queue.cur_time, 10);
    ASSERT_EQ(queue.to_vector(), evs({11}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, false);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{11});

    tracker.set_desired_event(ev(12), queue);
    ASSERT_EQ(queue.cur_time, 10);
    ASSERT_EQ(queue.to_vector(), evs({11}));
    ASSERT_EQ(tracker.has_queued_time, true);
    ASSERT_EQ(tracker.has_desired_time, true);
    ASSERT_EQ(tracker.queued_time, cyclic_time_int{11});
    ASSERT_EQ(tracker.desired_time, cyclic_time_int{12});
}

TEST(bit_bucket_queue, wraparound_all_the_way_around) {
    std::priority_queue<int64_t> reference_queue;
    bit_bucket_queue<true> actual_queue;
    for (size_t k = 0; k < 100; k++) {
        actual_queue.enqueue(TentativeEvent{cyclic_time_int{k}});
        reference_queue.push(-(int64_t)k);
    }

    std::mt19937 rng(0); // NOLINT(cert-msc51-cpp)
    size_t n = 0;
    while (!reference_queue.empty()) {
        auto actual = (size_t)-reference_queue.top();
        auto t = (size_t)actual_queue.dequeue().time.widen_from_nearby_reference(actual_queue.cur_time);
        ASSERT_EQ(t, actual) << n;
        reference_queue.pop();
        n++;
        if (n < (1 << 17)) {
            t += rng() % 1000;
            actual_queue.enqueue(TentativeEvent{cyclic_time_int{t}});
            reference_queue.push(-(int64_t)t);
        }
    }
}