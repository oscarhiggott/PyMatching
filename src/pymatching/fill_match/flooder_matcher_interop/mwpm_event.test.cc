#include "pymatching/fill_match/flooder_matcher_interop/mwpm_event.h"

#include <gtest/gtest.h>

#include "pymatching/fill_match/flooder/graph.h"
#include "pymatching/fill_match/flooder/graph_fill_region.h"
#include "pymatching/fixed_length_vector.h"

using namespace pm;

TEST(mwpm_event, basic_usage) {
    MwpmEvent ev{RegionHitRegionEventData{
        nullptr,
        nullptr,
        CompressedEdge{nullptr, nullptr, 2},
    }};
    ASSERT_EQ(
        ev.str(),
        "MwpmEvent{.type=REGION_HIT_REGION, .dat={.region1=0, .region2=0, .edge=CompressedEdge{.obs_mask=2, "
        ".loc_from=0, .loc_to=0}}}");
}
