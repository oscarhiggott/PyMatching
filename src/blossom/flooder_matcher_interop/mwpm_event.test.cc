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

#include "pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.h"

#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"

using namespace pm;

TEST(mwpm_event, basic_usage) {
    MwpmEvent ev{RegionHitRegionEventData{
        nullptr,
        nullptr,
        CompressedEdge{nullptr, nullptr, 2},
    }};
    ASSERT_TRUE(ev.str().starts_with("MwpmEvent{.type=REGION_HIT_REGION, .dat={.region1="));
}
