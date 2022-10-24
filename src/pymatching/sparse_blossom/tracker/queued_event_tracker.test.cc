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

#include "pymatching/sparse_blossom/tracker/queued_event_tracker.h"

#include <gtest/gtest.h>
#include <random>

using namespace pm;

TEST(QueuedEventTracker, basic_usage) {
    auto ev = [](int x) {
        return FloodCheckEvent{cyclic_time_int{x}};
    };
    auto evs = [&](std::vector<int> x) {
        std::vector<FloodCheckEvent> r;
        for (auto e : x) {
            r.push_back(ev(e));
        }
        return r;
    };
    radix_heap_queue<true> queue;
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
