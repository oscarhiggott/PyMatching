#include "gtest/gtest.h"
#include "stim.h"
#include "mwpm_decoding.h"
#include <fstream>


TEST(MwpmDecoding, CompareSolutionWeights) {
    auto shots_in = std::fopen("../data/surface_code_rotated_memory_x_13_0.01_1000_shots.b8", "r");
    auto dem_file = std::fopen("../data/surface_code_rotated_memory_x_13_0.01.dem", "r");

    ASSERT_TRUE(shots_in);
    ASSERT_TRUE(dem_file);
    stim::DetectorErrorModel dem = stim::DetectorErrorModel::from_file(dem_file);
    fclose(dem_file);

    size_t num_obs = dem.count_observables();
    std::unique_ptr<stim::MeasureRecordReader> obs_reader;
    auto reader = stim::MeasureRecordReader::make(shots_in,
                                                  stim::SAMPLE_FORMAT_B8,
                                                  0,
                                                  dem.count_detectors(),
                                                  dem.count_observables());

    pm::weight_int num_buckets = 1000;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    std::ifstream is("../data/surface_code_rotated_memory_x_13_0.01_1000_shots_1000_"
                     "buckets_solution_weights_pymatchingv0.7_exact.txt");
    std::istream_iterator<int> start(is), end;
    std::vector<int> expected_weights(start, end);

    std::ifstream is2("../data/surface_code_rotated_memory_x_13_0.01_1000_shots_1000_"
                      "buckets_predictions_pymatchingv0.7_exact.txt");
    std::istream_iterator<int> start2(is2), end2;
    std::vector<int> expected_obs_masks(start2, end2);

    stim::SparseShot sparse_shot;
    size_t num_mistakes = 0;
    size_t num_shots = 0;
    size_t max_shots = 500;
    while (reader->start_and_read_entire_record(sparse_shot)) {
        if (num_shots > max_shots)
            break;
        auto res = pm::decode_detection_events(mwpm, sparse_shot.hits);
        if (sparse_shot.obs_mask != res.obs_mask) {
            num_mistakes++;
        }
        ASSERT_EQ(res.weight, expected_weights[num_shots]);
        // Observable masks do not need to match exactly due to degeneracy, but they do for this dataset
        ASSERT_EQ(res.obs_mask, expected_obs_masks[num_shots]);
        sparse_shot.clear();
        num_shots++;
    }
}