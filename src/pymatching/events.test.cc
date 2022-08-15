#include "pymatching/events.h"

#include <gtest/gtest.h>

#include "pymatching/fixed_length_vector.h"

TEST(Events, TentativeEvent) {
    pm::TentativeEvent te{};
    te.tentative_event_type = pm::LOOK_AT_SHRINKING_REGION;
}
