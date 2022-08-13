#include "pymatching/mwpm_decoding.h"

#include "pymatching/perf/util.perf.h"
#include "stim.h"

std::pair<stim::DetectorErrorModel, std::vector<stim::SparseShot>> generate_data(
    size_t distance, size_t rounds, double noise, size_t num_shots) {
    stim::CircuitGenParameters gen(rounds, distance, "rotated_memory_x");
    gen.after_clifford_depolarization = noise;
    gen.after_reset_flip_probability = noise;
    gen.before_measure_flip_probability = noise;
    stim::Circuit circuit = stim::generate_surface_code_circuit(gen).circuit;
    std::mt19937_64 rng(0);  // NOLINT(cert-msc51-cpp)
    size_t num_detectors = circuit.count_detectors();
    auto results = stim::detector_samples(circuit, num_shots, false, true, rng);
    std::vector<stim::SparseShot> shots;
    for (size_t k = 0; k < num_shots; k++) {
        shots.push_back({});
        shots.back().obs_mask = results[num_detectors][k];
        for (size_t d = 0; d < num_detectors; d++) {
            if (results[d][k]) {
                shots.back().hits.push_back(d);
            }
        }
    }

    auto dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false);

    return {dem, shots};
}

BENCHMARK(Decode_surface_r11_d11_p100) {
    auto data = generate_data(11, 11, 0.01, 128);
    const auto &dem = data.first;
    const auto &shots = data.second;

    size_t num_buckets = 1024;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    size_t num_dets = 0;
    for (const auto &shot : shots) {
        num_dets += shot.hits.size();
    }

    size_t num_mistakes = 0;
    benchmark_go([&]() {
        for (const auto &shot : shots) {
            auto res = pm::decode_detection_events(mwpm, shot.hits);
            if (shot.obs_mask != res.obs_mask) {
                num_mistakes++;
            }
        }
    })
        .goal_millis(120)
        .show_rate("dets", (double)num_dets)
        .show_rate("shots", (double)shots.size());
    if (num_mistakes == shots.size()) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(Decode_surface_r11_d11_p1000) {
    auto data = generate_data(11, 11, 0.001, 512);
    const auto &dem = data.first;
    const auto &shots = data.second;

    size_t num_buckets = 1024;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    size_t num_dets = 0;
    for (const auto &shot : shots) {
        num_dets += shot.hits.size();
    }

    size_t num_mistakes = 0;
    benchmark_go([&]() {
        for (const auto &shot : shots) {
            auto res = pm::decode_detection_events(mwpm, shot.hits);
            if (shot.obs_mask != res.obs_mask) {
                num_mistakes++;
            }
        }
    })
        .goal_millis(22)
        .show_rate("dets", (double)num_dets)
        .show_rate("shots", (double)shots.size());
    if (num_mistakes == shots.size()) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(Decode_surface_r11_d11_p10000) {
    auto data = generate_data(11, 11, 0.0001, 256);
    const auto &dem = data.first;
    const auto &shots = data.second;

    size_t num_buckets = 1024;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    size_t num_dets = 0;
    for (const auto &shot : shots) {
        num_dets += shot.hits.size();
    }

    size_t num_mistakes = 0;
    benchmark_go([&]() {
        for (const auto &shot : shots) {
            auto res = pm::decode_detection_events(mwpm, shot.hits);
            if (shot.obs_mask != res.obs_mask) {
                num_mistakes++;
            }
        }
    })
        .goal_micros(800)
        .show_rate("dets", (double)num_dets)
        .show_rate("shots", (double)shots.size());
    if (num_mistakes == shots.size()) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(Decode_surface_r11_d11_p100000) {
    auto data = generate_data(11, 11, 0.00001, 1024);
    const auto &dem = data.first;
    const auto &shots = data.second;

    size_t num_buckets = 1024;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    size_t num_dets = 0;
    for (const auto &shot : shots) {
        num_dets += shot.hits.size();
    }

    size_t num_mistakes = 0;
    benchmark_go([&]() {
        for (const auto &shot : shots) {
            auto res = pm::decode_detection_events(mwpm, shot.hits);
            if (shot.obs_mask != res.obs_mask) {
                num_mistakes++;
            }
        }
    })
        .goal_micros(200)
        .show_rate("dets", (double)num_dets)
        .show_rate("shots", (double)shots.size());
    if (num_mistakes == shots.size()) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(Decode_surface_r21_d21_p100) {
    auto data = generate_data(21, 21, 0.01, 8);
    const auto &dem = data.first;
    const auto &shots = data.second;

    size_t num_buckets = 1024;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    size_t num_dets = 0;
    for (const auto &shot : shots) {
        num_dets += shot.hits.size();
    }

    size_t num_mistakes = 0;
    benchmark_go([&]() {
        for (const auto &shot : shots) {
            auto res = pm::decode_detection_events(mwpm, shot.hits);
            if (shot.obs_mask != res.obs_mask) {
                num_mistakes++;
            }
        }
    })
        .goal_millis(100)
        .show_rate("dets", (double)num_dets)
        .show_rate("shots", (double)shots.size());
    if (num_mistakes == 0) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(Decode_surface_r21_d21_p1000) {
    auto data = generate_data(21, 21, 0.001, 256);
    const auto &dem = data.first;
    const auto &shots = data.second;

    size_t num_buckets = 1024;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    size_t num_dets = 0;
    for (const auto &shot : shots) {
        num_dets += shot.hits.size();
    }

    size_t num_mistakes = 0;
    benchmark_go([&]() {
        for (const auto &shot : shots) {
            auto res = pm::decode_detection_events(mwpm, shot.hits);
            if (shot.obs_mask != res.obs_mask) {
                num_mistakes++;
            }
        }
    })
        .goal_millis(80)
        .show_rate("dets", (double)num_dets)
        .show_rate("shots", (double)shots.size());
    if (num_mistakes == shots.size()) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(Decode_surface_r21_d21_p10000) {
    auto data = generate_data(21, 21, 0.0001, 512);
    const auto &dem = data.first;
    const auto &shots = data.second;

    size_t num_buckets = 1024;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    size_t num_dets = 0;
    for (const auto &shot : shots) {
        num_dets += shot.hits.size();
    }

    size_t num_mistakes = 0;
    benchmark_go([&]() {
        for (const auto &shot : shots) {
            auto res = pm::decode_detection_events(mwpm, shot.hits);
            if (shot.obs_mask != res.obs_mask) {
                num_mistakes++;
            }
        }
    })
        .goal_millis(13)
        .show_rate("dets", (double)num_dets)
        .show_rate("shots", (double)shots.size());
    if (num_mistakes == shots.size()) {
        std::cerr << "data dependence";
    }
}

BENCHMARK(Decode_surface_r21_d21_p100000) {
    auto data = generate_data(21, 21, 0.00001, 512);
    const auto &dem = data.first;
    const auto &shots = data.second;

    size_t num_buckets = 1024;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    size_t num_dets = 0;
    for (const auto &shot : shots) {
        num_dets += shot.hits.size();
    }

    size_t num_mistakes = 0;
    benchmark_go([&]() {
        for (const auto &shot : shots) {
            auto res = pm::decode_detection_events(mwpm, shot.hits);
            if (shot.obs_mask != res.obs_mask) {
                num_mistakes++;
            }
        }
    })
        .goal_micros(900)
        .show_rate("dets", (double)num_dets)
        .show_rate("shots", (double)shots.size());
    if (num_mistakes == shots.size()) {
        std::cerr << "data dependence";
    }
}
