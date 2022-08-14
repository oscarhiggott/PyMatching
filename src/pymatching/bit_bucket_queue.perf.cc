#include "pymatching/bit_bucket_queue.h"
#include <random>
#include <iostream>

#include "pymatching/perf/util.perf.h"

BENCHMARK(bucket_queue_sort) {
    std::mt19937 rng(0); // NOLINT(cert-msc51-cpp)

    std::vector<pm::time_int> v;
    for (size_t k = 0; k < 1000; k++) {
        v.push_back((pm::time_int)(rng() & 0x3FFFFFFF));
    }

    bool dependence = false;
    benchmark_go([&]() {
        pm::bit_bucket_queue<false> q;
        for (pm::time_int t : v) {
            q.enqueue(pm::TentativeEvent(t));
        }
        pm::TentativeEvent out{};
        while (q.try_pop_valid(&out)) {
            q.force_pop_valid();
        }
        if (out.time == 0) {
            dependence = true;
        }
    }).goal_micros(60).show_rate("EnqueueDequeues", (double)v.size());
    if (dependence) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(bucket_queue_stream) {
    size_t n = 10000;

    bool dependence = false;
    benchmark_go([&]() {
        pm::bit_bucket_queue<false> q;
        for (size_t k = 0; k < 10; k++) {
            for (size_t r = 0; r < k; r++) {
                q.enqueue(pm::TentativeEvent((pm::time_int)k));
            }
        }
        for (size_t k = 0; k < n; k++) {
            q.enqueue(pm::TentativeEvent(q.force_pop_valid().time));
        }
    }).goal_micros(150).show_rate("EnqueueDequeues", (double)n);
    if (dependence) {
        std::cerr << "data dependence";
    }
}
