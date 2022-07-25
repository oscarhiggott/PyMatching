#ifndef PYMATCHING2_VARYING_PERF_H
#define PYMATCHING2_VARYING_PERF_H

#include <chrono>
#include <vector>
#include "varying.h"
#include <tuple>


namespace pm {

    enum FUNC_TYPE {DIRECT, CONDITIONAL};

    std::tuple<int32_t , double> time_get_distance(FUNC_TYPE type) {
        const int NUM_REPEATS = 10000000;
        std::vector<pm::Varying32> v;
        std::vector<int32_t> times;
        v.reserve(NUM_REPEATS);
        int32_t total_distance = 0;
        times.reserve(NUM_REPEATS);
        for (int i = 0; i < NUM_REPEATS; i++) {
            int32_t y_intercept = i - (NUM_REPEATS / 2);
            int32_t slope = i % 3;
            v.emplace_back((y_intercept << 2) + slope);
            times.push_back(i * 3);
        }
        auto start = std::chrono::steady_clock::now();
        if (type == DIRECT){
            for (int i = 0; i < NUM_REPEATS; i++) {
                total_distance += v[i].get_distance_at_time(times[i]);
            }
        } else if (type == CONDITIONAL) {
//            for (int i = 0; i < NUM_REPEATS; i++) {
//                total_distance +=  v[i].get_distance_at_time_conditional(times[i]);
//            }
            throw std::runtime_error("No longer implemented");
        } else {
            throw std::runtime_error("Didn't recognise type");
        }

        auto end = std::chrono::steady_clock::now();
        auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
        auto microseconds_per_run = (double)micros / (double)NUM_REPEATS;
        return std::make_tuple(total_distance, microseconds_per_run);
    }
}



#endif //PYMATCHING2_VARYING_PERF_H
