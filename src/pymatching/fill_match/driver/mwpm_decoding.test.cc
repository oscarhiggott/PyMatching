#include "pymatching/fill_match/driver/mwpm_decoding.h"

#include <fstream>

#include "gtest/gtest.h"

#include "stim.h"

std::string find_test_data_file(const char *name) {
    std::vector<std::string> directories_to_check = {
        "data/",
        "../data/",
        "../../data/",
    };
    for (const auto &d : directories_to_check) {
        std::string path = d + name;
        FILE *f = fopen((d + name).c_str(), "r");
        if (f != nullptr) {
            fclose(f);
            return path;
        }
    }
    throw std::invalid_argument("Failed to find test data file " + std::string(name));
}


struct DecodingTestCase{
    std::vector<int> expected_weights;
    std::vector<int> expected_obs_masks;
    stim::DetectorErrorModel detector_error_model;
    std::unique_ptr<stim::MeasureRecordReader> reader;
};

DecodingTestCase load_surface_code_d13_p100_test_case(){
    auto shots_in = std::fopen(find_test_data_file("surface_code_rotated_memory_x_13_0.01_1000_shots.b8").c_str(), "r");
    auto dem_file = std::fopen(find_test_data_file("surface_code_rotated_memory_x_13_0.01.dem").c_str(), "r");

    assert(shots_in);
    assert(dem_file);
    stim::DetectorErrorModel dem = stim::DetectorErrorModel::from_file(dem_file);
    fclose(dem_file);

    auto reader = stim::MeasureRecordReader::make(
            shots_in, stim::SAMPLE_FORMAT_B8, 0, dem.count_detectors(), dem.count_observables());


    std::ifstream is(
            find_test_data_file(
                    "surface_code_rotated_memory_x_13_0.01_1000_shots_1000_buckets_solution_weights_pymatchingv0.7_exact.txt")
                    .c_str());
    std::istream_iterator<int> start(is), end;
    std::vector<int> expected_weights(start, end);

    std::ifstream is2(
            find_test_data_file(
                    "surface_code_rotated_memory_x_13_0.01_1000_shots_1000_buckets_predictions_pymatchingv0.7_exact.txt")
                    .c_str());
    std::istream_iterator<int> start2(is2), end2;
    std::vector<int> expected_obs_masks(start2, end2);
    return {expected_weights, expected_obs_masks, dem, std::move(reader)};
}


TEST(MwpmDecoding, CompareSolutionWeights) {
    auto test_case = load_surface_code_d13_p100_test_case();

    pm::weight_int num_distinct_weights = 1001;
    auto mwpm = pm::detector_error_model_to_mwpm(test_case.detector_error_model, num_distinct_weights);

    stim::SparseShot sparse_shot;
    size_t num_mistakes = 0;
    size_t num_shots = 0;
    size_t max_shots = 500;
    while (test_case.reader->start_and_read_entire_record(sparse_shot)) {
        if (num_shots > max_shots)
            break;
        auto res = pm::decode_detection_events(mwpm, sparse_shot.hits);
        if (sparse_shot.obs_mask != res.obs_mask) {
            num_mistakes++;
        }
        EXPECT_EQ(res.weight, test_case.expected_weights[num_shots]);
        // Observable masks do not need to match exactly due to degeneracy, but they do for this dataset
        ASSERT_EQ(res.obs_mask, test_case.expected_obs_masks[num_shots]);
        sparse_shot.clear();
        num_shots++;
    }
}


pm::obs_int bit_vector_to_int_mask(std::vector<uint8_t> observables) {
    pm::obs_int obs = 0;
    for (size_t i = 0; i < observables.size(); i++){
        if (observables[i])
            obs ^= 1 << i;
    }
    return obs;
}


TEST(MwpmDecoding, CompareSolutionWeightsForMoreThan64Observables) {
    auto test_case = load_surface_code_d13_p100_test_case();

    pm::weight_int num_distinct_weights = 1001;
    test_case.detector_error_model.append_logical_observable_instruction(stim::DemTarget::observable_id(128));
    auto mwpm = pm::detector_error_model_to_mwpm(test_case.detector_error_model, num_distinct_weights);

    stim::SparseShot sparse_shot;
    size_t num_mistakes = 0;
    size_t num_shots = 0;
    size_t max_shots = 500;
    while (test_case.reader->start_and_read_entire_record(sparse_shot)) {
        if (num_shots > max_shots)
            break;
        auto res = pm::decode_detection_events_with_more_than_64_observables(mwpm, sparse_shot.hits);
        if (sparse_shot.obs_mask != bit_vector_to_int_mask(res.obs_crossed)) {
            num_mistakes++;
        }
        EXPECT_EQ(res.weight, test_case.expected_weights[num_shots]);
        // Observable masks do not need to match exactly due to degeneracy, but they do for this dataset
        ASSERT_EQ(bit_vector_to_int_mask(res.obs_crossed), test_case.expected_obs_masks[num_shots]);
        sparse_shot.clear();
        num_shots++;
    }
}
