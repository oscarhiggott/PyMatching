#include "pymatching/fill_match/driver/mwpm_decoding.h"

#include <fstream>

#include "gtest/gtest.h"

#include "pymatching/fill_match/flooder/graph.h"
#include "stim.h"

std::string find_test_data_file(const char* name) {
    std::vector<std::string> directories_to_check = {
        "data/",
        "../data/",
        "../../data/",
    };
    for (const auto& d : directories_to_check) {
        std::string path = d + name;
        FILE* f = fopen((d + name).c_str(), "r");
        if (f != nullptr) {
            fclose(f);
            return path;
        }
    }
    throw std::invalid_argument("Failed to find test data file " + std::string(name));
}

struct DecodingTestCase {
    std::vector<int> expected_weights;
    std::vector<int> expected_obs_masks;
    stim::DetectorErrorModel detector_error_model;
    std::unique_ptr<stim::MeasureRecordReader> reader;
};

DecodingTestCase load_surface_code_d13_p100_test_case() {
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
        auto res = pm::decode_detection_events_for_up_to_64_observables(mwpm, sparse_shot.hits);
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

TEST(MwpmDecoding, CompareSolutionWeightsWithNoLimitOnNumObservables) {
    for (size_t i : {0, 1}) {
        auto test_case = load_surface_code_d13_p100_test_case();
        pm::weight_int num_distinct_weights = 1001;
        if (i == 1)
            test_case.detector_error_model.append_logical_observable_instruction(stim::DemTarget::observable_id(128));
        auto mwpm = pm::detector_error_model_to_mwpm(test_case.detector_error_model, num_distinct_weights);

        stim::SparseShot sparse_shot;
        size_t num_mistakes = 0;
        size_t num_shots = 0;
        size_t max_shots = 100;
        pm::ExtendedMatchingResult res(mwpm.flooder.graph.num_observables);
        while (test_case.reader->start_and_read_entire_record(sparse_shot)) {
            if (num_shots > max_shots)
                break;
            pm::decode_detection_events(mwpm, sparse_shot.hits, res.obs_crossed.data(), res.weight);
            if (sparse_shot.obs_mask != res.obs_crossed[0]) {
                num_mistakes++;
            }
            EXPECT_EQ(res.weight, test_case.expected_weights[num_shots]);
            // Observable masks do not need to match exactly due to degeneracy, but they do for this dataset
            ASSERT_EQ(res.obs_crossed[0], test_case.expected_obs_masks[num_shots]);
            sparse_shot.clear();
            num_shots++;
            std::fill(res.obs_crossed.begin(), res.obs_crossed.end(), 0);
            res.weight = 0;
        }
    }
}

TEST(MwpmDecoding, DecodeToMatchEdges) {
    auto test_case = load_surface_code_d13_p100_test_case();

    pm::weight_int num_distinct_weights = 1001;
    auto mwpm = pm::detector_error_model_to_mwpm(test_case.detector_error_model, num_distinct_weights);

    stim::SparseShot sparse_shot;
    size_t num_shots = 0;
    size_t max_shots = 10;
    while (test_case.reader->start_and_read_entire_record(sparse_shot)) {
        if (num_shots > max_shots)
            break;
        pm::decode_detection_events_to_match_edges(mwpm, sparse_shot.hits);
        auto& match_edges = mwpm.flooder.match_edges;
        uint64_t obs_mask = 0;
        std::vector<uint64_t> dets;
        for (auto& e : match_edges) {
            obs_mask ^= e.obs_mask;
            dets.push_back(e.loc_from - &mwpm.flooder.graph.nodes[0]);
            if (e.loc_to)
                dets.push_back(e.loc_to - &mwpm.flooder.graph.nodes[0]);
        }

        // sparse_shot.hits already sorted. Compare with sorted dets.
        std::sort(dets.begin(), dets.end());
        ASSERT_EQ(dets, sparse_shot.hits);

        // Observable masks do not need to match exactly due to degeneracy, but they do for this dataset
        ASSERT_EQ(obs_mask, test_case.expected_obs_masks[num_shots]);
        sparse_shot.clear();
        num_shots++;
    }
}

TEST(MwpmDecoding, CompareSolutionObsWithMaxNumBuckets) {
    for (size_t i : {0, 1}) {
        auto test_case = load_surface_code_d13_p100_test_case();
        pm::weight_int num_distinct_weights = 1 << (sizeof(pm::weight_int) * 8 - 4);
        if (i == 1)
            test_case.detector_error_model.append_logical_observable_instruction(stim::DemTarget::observable_id(128));
        auto mwpm = pm::detector_error_model_to_mwpm(test_case.detector_error_model, num_distinct_weights);

        stim::SparseShot sparse_shot;
        size_t num_mistakes = 0;
        size_t num_shots = 0;
        size_t max_shots = 100;
        pm::ExtendedMatchingResult res(mwpm.flooder.graph.num_observables);
        while (test_case.reader->start_and_read_entire_record(sparse_shot)) {
            if (num_shots > max_shots)
                break;
            pm::decode_detection_events(mwpm, sparse_shot.hits, res.obs_crossed.data(), res.weight);
            if (sparse_shot.obs_mask != res.obs_crossed[0]) {
                num_mistakes++;
            }
            ASSERT_EQ(res.obs_crossed[0], test_case.expected_obs_masks[num_shots]);
            sparse_shot.clear();
            num_shots++;
            std::fill(res.obs_crossed.begin(), res.obs_crossed.end(), 0);
            res.weight = 0;
        }
    }
}

TEST(MwpmDecoding, FillBitVectorFromObsMask) {
    std::vector<uint8_t> expected_bit_vector = {0, 0, 0, 1, 0, 0, 0, 1, 0, 1};
    std::vector<uint8_t> bit_vector(10);
    pm::fill_bit_vector_from_obs_mask(648, bit_vector.data(), 10);
    ASSERT_EQ(bit_vector, expected_bit_vector);
}

TEST(MwpmDecoding, BitVectorToObsMask) {
    ASSERT_EQ(pm::bit_vector_to_obs_mask({0, 1, 0, 0, 0, 0, 1, 1, 0, 0}), 194);
}

TEST(MwpmDecoding, HandleAllNegativeWeights) {
    for (size_t num_nodes : {50, 80}) {
        auto mwpm = pm::Mwpm(
            pm::GraphFlooder(pm::MatchingGraph(num_nodes, num_nodes)), pm::SearchFlooder(pm::SearchGraph(num_nodes)));
        auto& g = mwpm.flooder.graph;
        for (size_t i = 0; i < num_nodes; i++)
            g.add_edge(i, (i + 1) % num_nodes, -2, {i});

        if (num_nodes > sizeof(pm::obs_int) * 8) {
            for (size_t i = 0; i < num_nodes; i++)
                mwpm.search_flooder.graph.add_edge(i, (i + 1) % num_nodes, -2, {i});
        }

        mwpm.flooder.sync_negative_weight_observables_and_detection_events();

        pm::ExtendedMatchingResult res(num_nodes);
        pm::decode_detection_events(mwpm, {10, 20}, res.obs_crossed.data(), res.weight);

        pm::ExtendedMatchingResult res_expected(num_nodes);
        for (size_t i = 0; i < num_nodes; i++) {
            if (i < 10 || i >= 20)
                res_expected.obs_crossed[i] ^= 1;
        }
        res_expected.weight = ((pm::signed_weight_int)num_nodes - 10) * -2;

        ASSERT_EQ(res, res_expected);

        if (num_nodes <= sizeof(pm::obs_int) * 8) {
            auto res2 = pm::decode_detection_events_for_up_to_64_observables(mwpm, {10, 20});
            ASSERT_EQ(res2.weight, res_expected.weight);
            pm::obs_int expected_obs_mask = 0;
            for (size_t i = 0; i < res_expected.obs_crossed.size(); i++) {
                if (res_expected.obs_crossed[i])
                    expected_obs_mask ^= (pm::obs_int)1 << i;
            }
            ASSERT_EQ(res2.obs_mask, expected_obs_mask);
        }
    }
}

TEST(MwpmDecoding, HandleSomeNegativeWeights) {
    size_t num_nodes = 8;
    for (size_t max_obs : {40, 100}) {
        auto mwpm = pm::Mwpm(
            pm::GraphFlooder(pm::MatchingGraph(num_nodes, max_obs + 1)), pm::SearchFlooder(pm::SearchGraph(num_nodes)));

        auto& g = mwpm.flooder.graph;
        g.add_boundary_edge(0, -4, {max_obs});
        for (size_t i = 0; i < 7; i += 2)
            g.add_edge(i, i + 1, 2, {i + 1});
        for (size_t i = 1; i < 7; i += 2)
            g.add_edge(i, i + 1, -4, {i + 1});
        g.add_boundary_edge(7, 2, {num_nodes});

        if (max_obs > sizeof(pm::obs_int) * 8) {
            auto& h = mwpm.search_flooder.graph;
            h.add_boundary_edge(0, -4, {max_obs});
            for (size_t i = 0; i < 7; i += 2)
                h.add_edge(i, i + 1, 2, {i + 1});
            for (size_t i = 1; i < 7; i += 2)
                h.add_edge(i, i + 1, -4, {i + 1});
            h.add_boundary_edge(7, 2, {num_nodes});
        }

        mwpm.flooder.sync_negative_weight_observables_and_detection_events();

        pm::ExtendedMatchingResult res(max_obs + 1);
        pm::decode_detection_events(mwpm, {0, 1, 2, 5, 6, 7}, res.obs_crossed.data(), res.weight);

        pm::ExtendedMatchingResult res_expected(max_obs + 1);
        res_expected.obs_crossed[max_obs] ^= 1;
        res_expected.obs_crossed[2] ^= 1;
        res_expected.obs_crossed[6] ^= 1;
        res_expected.obs_crossed[num_nodes] ^= 1;
        res_expected.weight = -10;

        ASSERT_EQ(res, res_expected);

        if (max_obs + 1 <= sizeof(pm::obs_int) * 8) {
            auto res2 = pm::decode_detection_events_for_up_to_64_observables(mwpm, {0, 1, 2, 5, 6, 7});
            ASSERT_EQ(res2.weight, res_expected.weight);
            pm::obs_int expected_obs_mask = 0;
            for (size_t i = 0; i < res_expected.obs_crossed.size(); i++) {
                if (res_expected.obs_crossed[i])
                    expected_obs_mask ^= (pm::obs_int)1 << i;
            }
            ASSERT_EQ(res2.obs_mask, expected_obs_mask);
        }
    }
}

TEST(MwpmDecoding, NoValidSolutionForLineGraph) {
    size_t num_nodes = 5;
    auto mwpm = pm::Mwpm(pm::GraphFlooder(pm::MatchingGraph(num_nodes, num_nodes)));
    auto& g = mwpm.flooder.graph;
    for (size_t i = 0; i < num_nodes; i++)
        g.add_edge(i, (i + 1) % num_nodes, 2, {i});
    pm::ExtendedMatchingResult res(num_nodes);
    EXPECT_THROW(pm::decode_detection_events(mwpm, {0, 2, 3}, res.obs_crossed.data(), res.weight);
                 , std::invalid_argument);
}

TEST(MwpmDecoding, InvalidSyndromeForToricCode) {
    auto shots_in = std::fopen(find_test_data_file("toric_code_unrotated_memory_x_5_0.005_1000.b8").c_str(), "r");
    auto dem_file = std::fopen(find_test_data_file("toric_code_unrotated_memory_x_5_0.005.dem").c_str(), "r");

    assert(shots_in);
    assert(dem_file);
    stim::DetectorErrorModel dem = stim::DetectorErrorModel::from_file(dem_file);
    fclose(dem_file);
    size_t num_distinct_weights = 1000;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_distinct_weights);
    auto reader = stim::MeasureRecordReader::make(
        shots_in, stim::SAMPLE_FORMAT_B8, 0, dem.count_detectors(), dem.count_observables());

    stim::SparseShot sparse_shot;
    size_t num_shots = 0;
    size_t max_shots = 10;
    while (reader->start_and_read_entire_record(sparse_shot)) {
        if (num_shots > max_shots)
            break;
        auto& detection_events = sparse_shot.hits;
        // Flip detector 0. detection_events should already be sorted.
        if (detection_events.empty() || detection_events[0] != 0) {
            detection_events.push_back(0);
        } else if (!detection_events.empty() && detection_events[0] == 0) {
            detection_events.erase(detection_events.begin());
        }
        EXPECT_THROW(pm::decode_detection_events_for_up_to_64_observables(mwpm, detection_events);
                     , std::invalid_argument);
        sparse_shot.clear();
        num_shots++;
    }
}
