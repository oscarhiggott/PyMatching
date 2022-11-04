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

#include "pymatching/sparse_blossom/diagram/animation_main.h"

#include <filesystem>

#include "pymatching/sparse_blossom/diagram/mwpm_diagram.h"
#include "pymatching/sparse_blossom/driver/io.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/matcher/mwpm.h"
#include "stim.h"

int pm::main_animation(int argc, const char **argv) {
    stim::check_for_unknown_arguments(
        {
            "--dets_in",
            "--dets_in_format",
            "--dets_in_includes_appended_observables",
            "--out_dir",
            "--dem_in",
            "--held_frames_per_event",
            "--held_frames_at_start",
            "--held_frames_at_end",
            "--max_growth_between_frames",
            "--pixels_per_unit_length",
            "--quiet",
        },
        {},
        "animate",
        argc,
        argv);

    FILE *dets_file = stim::find_open_file_argument("--dets_in", nullptr, "r", argc, argv);
    FILE *dem_file = stim::find_open_file_argument("--dem_in", nullptr, "r", argc, argv);
    stim::FileFormatData dets_in_format =
        stim::find_enum_argument("--dets_in_format", nullptr, stim::format_name_to_enum_map, argc, argv);
    bool dets_includes_obs = stim::find_bool_argument("--dets_in_includes_appended_observables", argc, argv);
    bool quiet = stim::find_bool_argument("--quiet", argc, argv);
    size_t held_frames_per_event = stim::find_int64_argument("--held_frames_per_event", -1, 0, 1000000, argc, argv);
    size_t held_frames_at_start = stim::find_int64_argument("--held_frames_at_start", -1, 0, 1000000, argc, argv);
    size_t held_frames_at_end = stim::find_int64_argument("--held_frames_at_end", -1, 0, 1000000, argc, argv);
    size_t max_growth_between_frames =
        stim::find_int64_argument("--max_growth_between_frames", -1, 0, 1000000, argc, argv);
    size_t pixels_per_unit_length = stim::find_int64_argument("--pixels_per_unit_length", -1, 4, 2048, argc, argv);
    std::string out_dir = stim::require_find_argument("--out_dir", argc, argv);
    if (!out_dir.ends_with('/') && !out_dir.ends_with('\\')) {
        out_dir.push_back('/');
    }

    auto dem = stim::DetectorErrorModel::from_file(dem_file);
    auto dets_reader = stim::MeasureRecordReader::make(
        dets_file, dets_in_format.id, 0, dem.count_detectors(), dem.count_observables() * dets_includes_obs);
    stim::SparseShot shot;
    if (!dets_reader->start_and_read_entire_record(shot)) {
        throw std::invalid_argument("No shot in file.");
    }
    fclose(dem_file);
    fclose(dets_file);

    auto coords = pm::pick_coords_for_drawing_from_dem(dem, pixels_per_unit_length);
    auto mwpm = pm::detector_error_model_to_mwpm(dem, 1000);

    if (!std::filesystem::exists(out_dir)) {
        std::filesystem::create_directory(out_dir);
    }
    pm::write_animated_decoding_svg_frames(
        mwpm,
        coords.first,
        coords.second,
        shot.hits,
        out_dir,
        !quiet,
        held_frames_per_event,
        held_frames_at_start,
        held_frames_at_end,
        max_growth_between_frames);

    return EXIT_SUCCESS;
}
