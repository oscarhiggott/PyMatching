#include "mwpm_decoding.h"

pm::Mwpm pm::detector_error_model_to_mwpm(
    const stim::DetectorErrorModel& detector_error_model, pm::weight_int num_buckets) {
    auto matching_graph = pm::detector_error_model_to_matching_graph(detector_error_model, num_buckets);
    pm::GraphFlooder flooder(matching_graph);
    return pm::Mwpm(flooder);
}

pm::MatchingResult pm::decode_detection_events(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    // Add detection events
    for (auto& detection : detection_events) {
        if (detection >= mwpm.flooder.graph.nodes.size())
            throw std::invalid_argument("Detection event index too large");
        mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[detection]);
    }

    while (true) {
        auto event = mwpm.flooder.next_event();
        if (event.event_type == pm::NO_EVENT)
            break;
        mwpm.process_event(event);
    }
    pm::MatchingResult res;
    for (auto& i : detection_events) {
        if (mwpm.flooder.graph.nodes[i].region_that_arrived)
            res += mwpm.shatter_blossom_and_extract_matches(mwpm.flooder.graph.nodes[i].top_region());
    }
    return res;
}
