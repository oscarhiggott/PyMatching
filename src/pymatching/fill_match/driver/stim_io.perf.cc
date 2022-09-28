#include "stim_io.h"

#include "mwpm_decoding.h"
#include "pymatching/perf/util.perf.h"

stim::DetectorErrorModel generate_dem(size_t distance, size_t rounds, double noise) {
    stim::CircuitGenParameters gen(rounds, distance, "rotated_memory_x");
    gen.after_clifford_depolarization = noise;
    gen.after_reset_flip_probability = noise;
    gen.before_measure_flip_probability = noise;
    stim::Circuit circuit = stim::generate_surface_code_circuit(gen).circuit;
    std::mt19937_64 rng(0);  // NOLINT(cert-msc51-cpp)
    auto dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false);
    return dem;
}

BENCHMARK(Load_dem_r11_d11_p100) {
    auto dem = generate_dem(11, 11, 0.01);
    size_t num_buckets = 1024;
    size_t num_loads = 10;
    benchmark_go([&]() {
        for (int i = 0; i < num_loads; i++)
            auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);
    })
        .goal_millis(35)
        .show_rate("loads", (double)num_loads);
}

BENCHMARK(Load_dem_r21_d21_p100) {
    auto dem = generate_dem(21, 21, 0.01);
    size_t num_buckets = 1024;
    size_t num_loads = 10;
    benchmark_go([&]() {
        for (int i = 0; i < num_loads; i++)
            auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);
    })
        .goal_millis(280)
        .show_rate("loads", (double)num_loads);
}
