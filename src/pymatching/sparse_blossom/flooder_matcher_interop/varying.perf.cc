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

#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"

#include <iostream>
#include <random>

#include "pymatching/perf/util.perf.h"

BENCHMARK(Varying32_get_distance_at_time) {
    const int NUM_ITEMS = 10000000;

    std::mt19937 rng(0);  // NOLINT(cert-msc51-cpp)
    std::vector<pm::Varying<int32_t>> varyings;
    std::vector<int32_t> times;
    varyings.resize(NUM_ITEMS);
    times.resize(NUM_ITEMS);
    for (size_t k = 0; k < NUM_ITEMS; k++) {
        varyings[k].data = (int32_t)(uint32_t)rng();
        times[k] = (int32_t)(uint32_t)rng();
    }
    size_t total = 0;
    benchmark_go([&]() {
        for (size_t k = 0; k < NUM_ITEMS; k++) {
            total += varyings[k].get_distance_at_time(times[k]);
        }
    })
        .goal_millis(3.5)
        .show_rate("calls", NUM_ITEMS);
    if (total == 0) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(Varying64_get_distance_at_time) {
    const int NUM_ITEMS = 10000000;

    std::mt19937 rng(0);  // NOLINT(cert-msc51-cpp)
    std::vector<pm::Varying<int64_t>> varyings;
    std::vector<int64_t> times;
    varyings.resize(NUM_ITEMS);
    times.resize(NUM_ITEMS);
    for (size_t k = 0; k < NUM_ITEMS; k++) {
        varyings[k].data = (int64_t)(uint64_t)rng();
        times[k] = (int64_t)(uint64_t)rng();
    }
    size_t total = 0;
    benchmark_go([&]() {
        for (size_t k = 0; k < NUM_ITEMS; k++) {
            total += varyings[k].get_distance_at_time(times[k]);
        }
    })
        .goal_millis(5.9)
        .show_rate("calls", NUM_ITEMS);
    if (total == 0) {
        std::cerr << "data dependence";
    }
}
