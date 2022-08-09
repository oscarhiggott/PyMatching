#include <iostream>
#include <vector>
#include "varying.h"
#include "stim.h"
#include "stim_io.h"


int main() {
    auto params = stim::CircuitGenParameters(5, 5, "rotated_memory_x");
    double p = 0.01;
    params.after_clifford_depolarization = p;
    params.before_round_data_depolarization = p;
    params.before_measure_flip_probability = p;
    params.after_reset_flip_probability = p;
    auto circuit = stim::generate_surface_code_circuit(params);
    auto dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(
            circuit.circuit, true, true, false,
            false, false,
            false
    );
    auto graph = pm::detector_error_model_to_matching_graph(dem);
    std::cout << graph.nodes.size() << std::endl;
    std::cout << dem.count_detectors() << std::endl;
}
