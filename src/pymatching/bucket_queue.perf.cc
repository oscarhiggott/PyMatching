#include "pymatching/bucket_queue.h"
#include <random>
#include <iostream>
#include <algorithm>
#include <immintrin.h>

#include "pymatching/perf/util.perf.h"

BENCHMARK(bucket_queue_sort) {
    std::mt19937 rng(0); // NOLINT(cert-msc51-cpp)

    std::vector<pm::time_int> v;
    for (size_t k = 0; k < 1000; k++) {
        v.push_back((pm::time_int)(rng() & 0x3FFFFFFF));
    }

    bool dependence = false;
    benchmark_go([&]() {
        pm::bit_bucket_queue q;
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
    }).goal_micros(70).show_rate("EnqueueDequeues", (double)v.size());
    if (dependence) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(bucket_queue_stream) {
    size_t n = 10000;

    bool dependence = false;
    benchmark_go([&]() {
        pm::bit_bucket_queue q;
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

BENCHMARK(compute_min_2) {
    uint16_t *values = (uint16_t *)_mm_malloc(sizeof(uint16_t) * 1024, 256);
    std::mt19937 rng(0);
    for (size_t k = 0; k < 1024; k++) {
        values[k] = rng();
    }
    size_t total = 0;
    benchmark_go([&]() {
        __m256i m = _mm256_set1_epi16(INT16_MIN);
        __m256i *p = (__m256i *)values;
        for (size_t k = 0; k < 1024 / 16; k++) {
            m = _mm256_min_epu16(p[k], m);
        }
        m = _mm256_min_epu16(m, _mm256_bsrli_epi128(m, 2));
        m = _mm256_min_epu16(m, _mm256_bsrli_epi128(m, 4));
        m = _mm256_min_epu16(m, _mm256_bsrli_epi128(m, 8));
        uint16_t m16 = _mm256_extract_epi16(m, 0);
        m16 = std::min(m16, (uint16_t)_mm256_extract_epi16(m, 8));
        total |= m16;
    });
    if (total == 0) {
        std::cerr << "data dependence";
    }
    std::cerr << total << "\n";
    delete values;
}
