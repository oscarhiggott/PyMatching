#include "gtest/gtest.h"
#include "graph_flooder.h"
#include "events.h"

TEST(GraphFlooder, PriorityQueue){
    pm::Graph graph(10);
    pm::GraphFillRegion gfr;
    pm::GraphFlooder flooder(graph);
    flooder.queue.emplace(&graph.nodes[0], 0, &graph.nodes[1], 1, 10);
    flooder.queue.emplace(&graph.nodes[1], 0, &graph.nodes[2], 1, 8);
    flooder.queue.emplace(&graph.nodes[2], 0, &graph.nodes[3], 1, 5);
    flooder.queue.emplace(&gfr, 70);
    flooder.queue.emplace(&graph.nodes[4], 0, &graph.nodes[5], 1, 100);
    auto e = flooder.queue.top();
    ASSERT_EQ(e.time, 5);
    ASSERT_EQ(e.tentative_event_type, pm::INTERACTION);
    ASSERT_EQ(e.neighbor_interaction_event_data.detector_node_1, &graph.nodes[2]);
    flooder.queue.pop();
    ASSERT_EQ(flooder.queue.top().time, 8);
    flooder.queue.pop();
    ASSERT_EQ(flooder.queue.top().time, 10);
    flooder.queue.pop();
    auto e2 = flooder.queue.top();
    ASSERT_EQ(e2.time, 70);
    ASSERT_EQ(e2.tentative_event_type, pm::SHRINKING);
    ASSERT_EQ(e2.region_shrink_event_data.region, &gfr);
    ASSERT_EQ(e2.is_invalidated, false);
    flooder.queue.pop();
    ASSERT_EQ(flooder.queue.top().time, 100);
}

