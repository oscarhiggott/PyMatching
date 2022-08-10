#include <iostream>
#include <vector>

#include "stim.h"
#include "stim_io.h"
#include "stim/simulators/detection_simulator.h"
#include "mwpm_decoding.h"
#include <chrono>
#include <cstring>


int main(int argc, char * argv[]) {
    if (argc < 4){
        std::cout << "Too few arguments" << std::endl;
        return 1;
    }
    auto f_dem = std::fopen(argv[1], "r");
    if (!f_dem){
        std::cout << "File " << argv[1] << " could not be opened" << std::endl;
        return 1;
    }
    auto f_shots = std::fopen(argv[2], "r");
    if (!f_shots){
        std::cout << "File " << argv[2] << " could not be opened" << std::endl;
        return 1;
    }

    auto dem = stim::DetectorErrorModel::from_file(f_dem);
    auto reader = stim::MeasureRecordReader::make(f_shots, /*stim::format_name_to_enum_map*/
                                                  stim::SAMPLE_FORMAT_B8,
                                                  0,
                                                  dem.count_detectors(),
                                                  dem.count_observables());
    pm::weight_int num_buckets = 10000;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    stim::SparseShot sparse_shot;
    sparse_shot.clear();
    size_t num_read = 0;
    size_t num_failures = 0;
    auto start = std::chrono::steady_clock::now();
    while (reader->start_and_read_entire_record(sparse_shot)) {
        auto res = pm::decode_detection_events(mwpm, sparse_shot.hits);
        if (std::strcmp(argv[3], "predictions") == 0) {
            std::cout << res.obs_mask << std::endl;
        } else if (std::strcmp(argv[3], "weight") == 0) {
            std::cout << res.weight << std::endl;
        }
        num_failures += res.obs_mask != sparse_shot.obs_mask;

        sparse_shot.clear();
        num_read++;
    }
    auto end = std::chrono::steady_clock::now();
    auto micros = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    auto microseconds_per_run = (double)micros / (double)num_read;
    if (std::strcmp(argv[3], "time") == 0)
        std::cout << microseconds_per_run << std::endl;
    if (std::strcmp(argv[3], "failures") == 0)
        std::cout << num_failures << " " << num_read << std::endl;
    std::fclose(f_dem);
    std::fclose(f_shots);
}
