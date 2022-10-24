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

#include "pymatching/sparse_blossom/driver/namespaced_main.h"

#include <chrono>
#include <cstring>
#include <iostream>
#include <vector>

#include "pymatching/sparse_blossom/diagram/animation_main.h"
#include "pymatching/sparse_blossom/driver/io.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "stim.h"
#include "stim/simulators/detection_simulator.h"

int main_predict(int argc, const char **argv) {
    stim::check_for_unknown_arguments(
        {
            "--in",
            "--in_format",
            "--in_includes_appended_observables",
            "--out",
            "--out_format",
            "--dem",
        },
        {},
        "predict",
        argc,
        argv);

    FILE *shots_in = stim::find_open_file_argument("--in", stdin, "r", argc, argv);
    FILE *predictions_out = stim::find_open_file_argument("--out", stdout, "w", argc, argv);
    FILE *dem_file = stim::find_open_file_argument("--dem", nullptr, "r", argc, argv);
    stim::FileFormatData shots_in_format =
        stim::find_enum_argument("--in_format", "b8", stim::format_name_to_enum_map, argc, argv);
    stim::FileFormatData predictions_out_format =
        stim::find_enum_argument("--out_format", "01", stim::format_name_to_enum_map, argc, argv);
    bool append_obs = stim::find_bool_argument("--in_includes_appended_observables", argc, argv);

    stim::DetectorErrorModel dem = stim::DetectorErrorModel::from_file(dem_file);
    fclose(dem_file);

    size_t num_obs = dem.count_observables();
    auto reader =
        stim::MeasureRecordReader::make(shots_in, shots_in_format.id, 0, dem.count_detectors(), append_obs * num_obs);
    auto writer = stim::MeasureRecordWriter::make(predictions_out, predictions_out_format.id);
    writer->begin_result_type('L');

    pm::weight_int num_buckets = pm::NUM_DISTINCT_WEIGHTS;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    stim::SparseShot sparse_shot;
    sparse_shot.clear();
    pm::ExtendedMatchingResult res(mwpm.flooder.graph.num_observables);
    while (reader->start_and_read_entire_record(sparse_shot)) {
        pm::decode_detection_events(mwpm, sparse_shot.hits, res.obs_crossed.data(), res.weight);
        for (size_t k = 0; k < num_obs; k++) {
            writer->write_bit(res.obs_crossed[k]);
        }
        writer->write_end();
        sparse_shot.clear();
        res.reset();
    }
    if (predictions_out != stdout) {
        fclose(predictions_out);
    }
    if (shots_in != stdin) {
        fclose(shots_in);
    }

    return EXIT_SUCCESS;
}

int main_count_mistakes(int argc, const char **argv) {
    stim::check_for_unknown_arguments(
        {
            "--in",
            "--in_format",
            "--in_includes_appended_observables",
            "--obs_in",
            "--obs_in_format",
            "--out",
            "--dem",
            "--time",
        },
        {},
        "count_mistakes",
        argc,
        argv);

    FILE *shots_in = stim::find_open_file_argument("--in", stdin, "r", argc, argv);
    FILE *obs_in = stim::find_open_file_argument("--obs_in", stdin, "r", argc, argv);
    FILE *stats_out = stim::find_open_file_argument("--out", stdout, "w", argc, argv);
    FILE *dem_file = stim::find_open_file_argument("--dem", nullptr, "r", argc, argv);
    stim::FileFormatData shots_in_format =
        stim::find_enum_argument("--in_format", "01", stim::format_name_to_enum_map, argc, argv);
    stim::FileFormatData obs_in_format =
        stim::find_enum_argument("--obs_in_format", "01", stim::format_name_to_enum_map, argc, argv);
    bool append_obs = stim::find_bool_argument("--in_includes_appended_observables", argc, argv);
    bool time = stim::find_bool_argument("--time", argc, argv);
    if (!append_obs && obs_in == nullptr) {
        throw std::invalid_argument("Must specify --in_includes_appended_observables or --obs_in.");
    }

    stim::DetectorErrorModel dem = stim::DetectorErrorModel::from_file(dem_file);
    fclose(dem_file);

    size_t num_obs = dem.count_observables();
    std::unique_ptr<stim::MeasureRecordReader> obs_reader;
    if (obs_in != stdin) {
        obs_reader = stim::MeasureRecordReader::make(obs_in, obs_in_format.id, 0, 0, num_obs);
    }
    auto reader =
        stim::MeasureRecordReader::make(shots_in, shots_in_format.id, 0, dem.count_detectors(), append_obs * num_obs);

    pm::weight_int num_buckets = pm::NUM_DISTINCT_WEIGHTS;
    auto mwpm = pm::detector_error_model_to_mwpm(dem, num_buckets);

    stim::SparseShot sparse_shot;
    stim::SparseShot obs_shot;
    size_t num_mistakes = 0;
    size_t num_shots = 0;
    auto start = std::chrono::steady_clock::now();
    while (reader->start_and_read_entire_record(sparse_shot)) {
        if (obs_reader == nullptr) {
            obs_shot.obs_mask = sparse_shot.obs_mask;
        } else {
            if (!obs_reader->start_and_read_entire_record(obs_shot)) {
                throw std::invalid_argument("Obs data ended before shot data ended.");
            }
        }
        auto res = pm::decode_detection_events_for_up_to_64_observables(mwpm, sparse_shot.hits);
        if (obs_shot.obs_mask != res.obs_mask) {
            num_mistakes++;
        }
        sparse_shot.clear();
        obs_shot.clear();
        num_shots++;
    }
    fprintf(stats_out, "%zu / %zu\n", num_mistakes, num_shots);
    if (stats_out != stdout) {
        fclose(stats_out);
    }
    if (shots_in != stdin) {
        fclose(shots_in);
    }

    auto end = std::chrono::steady_clock::now();
    auto microseconds = (double)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    if (time) {
        std::cerr << "Total decoding time: " << (int)microseconds << "us\n";
        std::cerr << "Decoding time per shot: " << (microseconds / num_shots) << "us\n";
    }

    return EXIT_SUCCESS;
}

int pm::main(int argc, const char **argv) {
    const char *command = "";
    if (argc >= 2) {
        command = argv[1];
    }
    try {
        if (strcmp(command, "predict") == 0) {
            return main_predict(argc, argv);
        }
        if (strcmp(command, "count_mistakes") == 0) {
            return main_count_mistakes(argc, argv);
        }
        if (strcmp(command, "animate") == 0) {
            return pm::main_animation(argc, argv);
        }
    } catch (std::invalid_argument &ex) {
        std::cerr << ex.what() << "\n";
        return EXIT_FAILURE;
    }

    std::stringstream ss;
    ss << "Unrecognized command. Available commands are:\n";
    ss << "    pymatching predict --dem file [--in file] [--out file] [--in_format 01|b8|...] [--out_format 01|b8|...] "
          "[--in_includes_appended_observables]\n";
    ss << "    pymatching count_mistakes --dem file [--in file] [--out file] [--in_format 01|b8|...] [--out_format "
          "01|B8|...] [--in_includes_appended_observables] [--obs_in] [--obs_in_format]\n";
    ss << "    pymatching animate "
          "--dets_in <file> "
          "--dets_in_format 01|b8|... "
          "--out_dir <directory> "
          "--dem_in <file> "
          "--held_frames_per_event # "
          "--held_frames_at_start # "
          "--held_frames_at_end # "
          "--max_growth_between_frames # "
          "--max_edge_weight # "
          "--pixels_per_unit_length # "
          "[--dets_in_includes_appended_observables] "
          "[--quiet]";
    throw std::invalid_argument(ss.str());
}
