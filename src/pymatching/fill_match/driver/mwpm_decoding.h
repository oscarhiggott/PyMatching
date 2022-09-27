#ifndef PYMATCHING2_MWPM_DECODING_H
#define PYMATCHING2_MWPM_DECODING_H

#include "pymatching/fill_match/driver/stim_io.h"
#include "pymatching/fill_match/matcher/mwpm.h"
#include "stim.h"

namespace pm {

std::vector<uint8_t> obs_mask_to_bit_vector(pm::obs_int obs_mask, size_t num_observables);

obs_int bit_vector_to_obs_mask(const std::vector<uint8_t>& bit_vector);

Mwpm detector_error_model_to_mwpm(
    const stim::DetectorErrorModel& detector_error_model, pm::weight_int num_distinct_weights);

pm::MatchingResult decode_detection_events(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

pm::ExtendedMatchingResult decode_detection_events_with_no_limit_on_num_observables(
        pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events
        );

}  // namespace pm

#endif  // PYMATCHING2_MWPM_DECODING_H
