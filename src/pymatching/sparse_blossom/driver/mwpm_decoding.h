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

#ifndef PYMATCHING2_MWPM_DECODING_H
#define PYMATCHING2_MWPM_DECODING_H

#include "pymatching/sparse_blossom/driver/io.h"
#include "pymatching/sparse_blossom/matcher/mwpm.h"
#include "stim.h"

namespace pm {

struct ExtendedMatchingResult {
    std::vector<uint8_t> obs_crossed;
    total_weight_int weight;
    ExtendedMatchingResult();
    explicit ExtendedMatchingResult(size_t num_observables);

    bool operator==(const ExtendedMatchingResult& rhs) const;

    bool operator!=(const ExtendedMatchingResult& rhs) const;

    ExtendedMatchingResult(std::vector<uint8_t> obs_crossed, total_weight_int weight);

    void reset();

    ExtendedMatchingResult& operator+=(const ExtendedMatchingResult& rhs);
    ExtendedMatchingResult operator+(const ExtendedMatchingResult& rhs) const;
};

inline void ExtendedMatchingResult::reset() {
    std::fill(obs_crossed.begin(), obs_crossed.end(), 0);
    weight = 0;
}

void fill_bit_vector_from_obs_mask(pm::obs_int obs_mask, uint8_t* obs_begin_ptr, size_t num_observables);
obs_int bit_vector_to_obs_mask(const std::vector<uint8_t>& bit_vector);

Mwpm detector_error_model_to_mwpm(
    const stim::DetectorErrorModel& detector_error_model,
    pm::weight_int num_distinct_weights,
    bool ensure_search_flooder_included = false);

MatchingResult decode_detection_events_for_up_to_64_observables(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

/// Used to decode detection events for an existing Mwpm object `mwpm', and a vector of
/// detection event indices `detection_events'. The predicted observables are XOR-ed into an
/// existing uint8_t array with at least `mwpm.flooder.graph.num_observables' elements,
/// the pointer to the first element of which is passed as the `obs_begin_ptr' argument.
/// The weight of the MWPM solution is added to the `weight' argument.
void decode_detection_events(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    uint8_t* obs_begin_ptr,
    pm::total_weight_int& weight);

/// Decode detection events using a Mwpm object and vector of detection event indices
/// Returns the compressed edges in the matching: the pairs of detection events that are
/// matched to each other via paths.
void decode_detection_events_to_match_edges(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

}  // namespace pm

#endif  // PYMATCHING2_MWPM_DECODING_H
