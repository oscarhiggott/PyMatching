#ifndef PYMATCHING2_MWPM_DECODING_H
#define PYMATCHING2_MWPM_DECODING_H

#include "pymatching/fill_match/driver/stim_io.h"
#include "pymatching/fill_match/matcher/mwpm.h"
#include "stim.h"

namespace pm {

void fill_bit_vector_from_obs_mask(pm::obs_int obs_mask, std::vector<uint8_t>::iterator obs_it_begin,
        std::vector<uint8_t>::iterator obs_it_end);
obs_int bit_vector_to_obs_mask(const std::vector<uint8_t>& bit_vector);

Mwpm detector_error_model_to_mwpm(
    const stim::DetectorErrorModel& detector_error_model, pm::weight_int num_distinct_weights);

MatchingResult decode_detection_events(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

void decode_detection_events_with_no_limit_on_num_observables(
        pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events,
        std::vector<uint8_t>::iterator obs_it_begin, pm::cumulative_time_int& weight);

}  // namespace pm

#endif  // PYMATCHING2_MWPM_DECODING_H
