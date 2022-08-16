#include "pymatching/fill_match/flooder_matcher_interop/mwpm_event.h"

#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>

#include "pymatching/fill_match/flooder/graph.h"
#include "pymatching/fill_match/flooder/graph_fill_region.h"

using namespace pm;
using ::testing::MatchesRegex;

TEST(mwpm_event, basic_usage) {
    MwpmEvent ev{RegionHitRegionEventData{
        nullptr,
        nullptr,
        CompressedEdge{nullptr, nullptr, 2},
    }};
    EXPECT_THAT(
            ev.str(),
            MatchesRegex("MwpmEvent{\\.type=REGION_HIT_REGION, \\.dat={\\.region1=0?x?0, \\.region2=0?x?0, "
                         "\\.edge=CompressedEdge{\\.obs_mask=2, \\.loc_from=0?x?0, \\.loc_to=0?x?0}}}")
            );
}
