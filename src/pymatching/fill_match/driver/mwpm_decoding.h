#ifndef PYMATCHING2_MWPM_DECODING_H
#define PYMATCHING2_MWPM_DECODING_H

#include "pymatching/fill_match/driver/stim_io.h"
#include "pymatching/fill_match/matcher/mwpm.h"
#include "stim.h"

namespace pm {

Mwpm detector_error_model_to_mwpm(
    const stim::DetectorErrorModel& detector_error_model, pm::weight_int num_distinct_weights);

pm::MatchingResult decode_detection_events(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

}  // namespace pm

#endif  // PYMATCHING2_MWPM_DECODING_H
