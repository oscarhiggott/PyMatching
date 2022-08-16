#include <gtest/gtest.h>

#include "pymatching/fill_match/flooder.h"

using namespace pm;

TEST(flood_check_event, basic_usage) {
    FloodCheckEvent ev(cyclic_time_int{6});

    ASSERT_EQ(ev.tentative_event_type, NO_TENTATIVE_EVENT);
    ASSERT_EQ(ev.time, 6);
    ASSERT_EQ(ev.str(), "TentativeEvent{.time=6, .type=NO_TENTATIVE_EVENT}");

    DetectorNode node;
    ev = {
        TentativeEventData_LookAtNode{&node},
        cyclic_time_int{5},
    };
    ASSERT_EQ(ev.tentative_event_type, LOOK_AT_NODE);
    ASSERT_EQ(ev.time, 5);
    ASSERT_EQ(ev.data_look_at_node.detector_node, &node);
    ASSERT_EQ(ev.str().find("TentativeEvent{.time=5, .type=LOOK_AT_NODE, .node="), 0);

    GraphFillRegion region;
    ev = {
        TentativeEventData_LookAtShrinkingRegion{&region},
        cyclic_time_int{4},
    };
    ASSERT_EQ(ev.tentative_event_type, LOOK_AT_SHRINKING_REGION);
    ASSERT_EQ(ev.time, 4);
    ASSERT_EQ(ev.data_look_at_shrinking_region.region, &region);
    ASSERT_EQ(ev.str().find("TentativeEvent{.time=4, .type=LOOK_AT_SHRINKING_REGION, .region="), 0);
}
