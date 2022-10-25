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

#include "pymatching/sparse_blossom/tracker/flood_check_event.h"

#include <gtest/gtest.h>

#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"

using namespace pm;

TEST(flood_check_event, basic_usage) {
    FloodCheckEvent ev(cyclic_time_int{6});

    ASSERT_EQ(ev.tentative_event_type, NO_FLOOD_CHECK_EVENT);
    ASSERT_EQ(ev.time, 6);
    ASSERT_EQ(ev.str(), "TentativeEvent{.time=6, .type=NO_TENTATIVE_EVENT}");

    DetectorNode node;
    ev = {
        &node,
        cyclic_time_int{5},
    };
    ASSERT_EQ(ev.tentative_event_type, LOOK_AT_NODE);
    ASSERT_EQ(ev.time, 5);
    ASSERT_EQ(ev.data_look_at_node, &node);
    ASSERT_EQ(ev.str().find("TentativeEvent{.time=5, .type=LOOK_AT_NODE, .node="), 0);

    GraphFillRegion region;
    ev = {
        &region,
        cyclic_time_int{4},
    };
    ASSERT_EQ(ev.tentative_event_type, LOOK_AT_SHRINKING_REGION);
    ASSERT_EQ(ev.time, 4);
    ASSERT_EQ(ev.data_look_at_shrinking_region, &region);
    ASSERT_EQ(ev.str().find("TentativeEvent{.time=4, .type=LOOK_AT_SHRINKING_REGION, .region="), 0);
}
