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

#include "pymatching/sparse_blossom/tracker/radix_heap_queue.h"

#include <iostream>
#include <random>

#include "pymatching/perf/util.perf.h"

using namespace pm;

BENCHMARK(bucket_queue_sort) {
    std::mt19937 rng(0);  // NOLINT(cert-msc51-cpp)

    std::vector<cyclic_time_int> v;
    for (size_t k = 0; k < 1000; k++) {
        v.push_back((cyclic_time_int)(rng() & 0x7FFF));
    }

    bool dependence = false;
    benchmark_go([&]() {
        radix_heap_queue<false> q;
        for (auto t : v) {
            q.enqueue(FloodCheckEvent(t));
        }
        while (true) {
            FloodCheckEvent out = q.dequeue();
            if (out.tentative_event_type == NO_FLOOD_CHECK_EVENT) {
                break;
            }
            if (out.time == 0) {
                dependence = true;
            }
        }
    })
        .goal_micros(4.9)
        .show_rate("EnqueueDequeues", (double)v.size());
    if (dependence) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(bucket_queue_stream) {
    size_t n = 10000;

    bool dependence = false;
    benchmark_go([&]() {
        radix_heap_queue<false> q;
        for (size_t k = 0; k < 10; k++) {
            for (size_t r = 0; r < k; r++) {
                q.enqueue(FloodCheckEvent(cyclic_time_int{k}));
            }
        }
        for (size_t k = 0; k < n; k++) {
            q.enqueue(FloodCheckEvent(q.dequeue().time + cyclic_time_int{100}));
        }
    })
        .goal_micros(99)
        .show_rate("EnqueueDequeues", (double)n);
    if (dependence) {
        std::cerr << "data dependence";
    }
}
