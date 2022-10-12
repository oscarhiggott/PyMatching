#include "pymatching/fill_match/flooder_matcher_interop/mwpm_event.h"

#include <gtest/gtest.h>

#include "pymatching/fill_match/flooder/graph.h"
#include "pymatching/fill_match/flooder/graph_fill_region.h"

using namespace pm;

TEST(mwpm_event, basic_usage) {
    MwpmEvent ev{RegionHitRegionEventData{
        nullptr,
        nullptr,
        CompressedEdge{nullptr, nullptr, 2},
    }};
    ASSERT_TRUE(ev.str().starts_with("MwpmEvent{.type=REGION_HIT_REGION, .dat={.region1="));
}
