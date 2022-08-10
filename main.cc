#include <iostream>
#include <vector>
#include "varying.h"
#include "stim.h"
#include "stim_io.h"
#include "stim/simulators/detection_simulator.h"
#include "mwpm_decoding.h"
#include <random>

#include <filesystem>
namespace fs = std::filesystem;


//std::mt19937_64 & global_urng() {
//    static std::mt19937_64 u{}; return u;
//}


int main() {
    auto f_circuit = std::fopen("../circuits/surface_code_rotated_memory_x_11_0.01.stim", "r");
    auto circuit = stim::Circuit::from_file(f_circuit);
    auto f_dem = std::fopen("../circuits/surface_code_rotated_memory_x_11_0.01.dem", "r");
    auto dem = stim::DetectorErrorModel::from_file(f_dem);
    auto f_shots = std::fopen("../circuits/surface_code_rotated_memory_x_11_0.01_50000.b8", "r");
    if (!f_shots)
        std::cout << "File not loaded" << std::endl;
    auto reader = stim::MeasureRecordReader::make(f_shots, /*stim::format_name_to_enum_map*/
                                                  stim::SAMPLE_FORMAT_B8,
                                                  0,
                                                  circuit.count_detectors(),
                                                  circuit.count_observables());
    pm::weight_int num_buckets = 1000;
    std::cout << dem.count_detectors() << " " << dem.count_detectors() << " " << std::endl;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    stim::SparseShot sparse_shot;
    sparse_shot.clear();
    size_t num_read = 0;
    size_t num_shots = 1000;
    while (num_read < num_shots && reader->start_and_read_entire_record(sparse_shot)) {
        auto res = pm::decode_detection_events(mwpm, sparse_shot.hits);
        std::cout << res.obs_mask << " " << sparse_shot.obs_mask << std::endl;
        sparse_shot.clear();
        num_read++;
    }
}
