#include <gtest/gtest.h>
#include "events.h"

TEST(Events, TentativeEvent) {
    pm::TentativeEvent te{};
    te.tentative_event_type = pm::SHRINKING;
}
