#include "pymatching/bucket_queue.h"
#include <random>
#include <iostream>
#include <algorithm>
#include <immintrin.h>

#include "pymatching/perf/util.perf.h"

BENCHMARK(compute_min_1) {
    uint16_t *values = (uint16_t *)_mm_malloc(sizeof(uint16_t) * 1024, 256);
    std::mt19937 rng(0);
    for (size_t k = 0; k < 1024; k++) {
        values[k] = rng();
    }
    size_t total = 0;
    benchmark_go([&]() {
        uint16_t m = UINT16_MAX;
        for (size_t k = 0; k < 1024; k++) {
            m = std::min(m, values[k]);
        }
        total |= m;
    });
    if (total == 0) {
        std::cerr << "data dependence";
    }
    std::cerr << total << "\n";
    delete values;
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
