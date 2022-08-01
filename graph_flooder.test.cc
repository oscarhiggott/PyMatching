#include "gtest/gtest.h"
#include "graph_flooder.h"
#include "events.h"

TEST(GraphFlooder, PriorityQueue){
    pm::Graph graph(10);
    pm::GraphFillRegion gfr;
    pm::GraphFlooder flooder(graph);
    auto e1 = new pm::TentativeEvent(&graph.nodes[0], 0, &graph.nodes[1], 1, 10);
    flooder.queue.push(e1);
    auto e2 = new pm::TentativeEvent(&graph.nodes[1], 0, &graph.nodes[2], 1, 8);
    flooder.queue.push(e2);
    auto e3 = new pm::TentativeEvent(&graph.nodes[2], 0, &graph.nodes[3], 1, 5);
    flooder.queue.push(e3);
    auto e4 = new pm::TentativeEvent(&gfr, 70);
    flooder.queue.push(e4);
    auto e5 = new pm::TentativeEvent(&graph.nodes[4], 0, &graph.nodes[5], 1, 100);
    flooder.queue.push(e5);
    auto e = flooder.queue.top();
    ASSERT_EQ(e->time, 5);
    ASSERT_EQ(e->tentative_event_type, pm::INTERACTION);
    ASSERT_EQ(e->neighbor_interaction_event_data.detector_node_1, &graph.nodes[2]);
    flooder.queue.pop();
    ASSERT_EQ(flooder.queue.top()->time, 8);
    flooder.queue.pop();
    ASSERT_EQ(flooder.queue.top()->time, 10);
    flooder.queue.pop();
    auto e6 = flooder.queue.top();
    ASSERT_EQ(e6->time, 70);
    ASSERT_EQ(e6->tentative_event_type, pm::SHRINKING);
    ASSERT_EQ(e6->region_shrink_event_data.region, &gfr);
    ASSERT_EQ(e6->is_invalidated, false);
    flooder.queue.pop();
    ASSERT_EQ(flooder.queue.top()->time, 100);
}


TEST(GraphFlooder, CreateRegion) {
    auto g = pm::Graph(5);
    g.add_boundary_edge(0, 3, 0);
    g.add_edge(0, 1, 5, 0);
    g.add_edge(1, 2, 11, 0);
    g.add_edge(2, 3, 100, 0);
    g.add_boundary_edge(3, 1000, 0);
    pm::GraphFlooder flooder(g);
    flooder.create_region(&flooder.graph.nodes[0]);
    flooder.create_region(&flooder.graph.nodes[2]);
    flooder.create_region(&flooder.graph.nodes[3]);
    auto e1 = flooder.queue.top();
    ASSERT_EQ(*e1, pm::TentativeEvent(
            &flooder.graph.nodes[0], 0,
            nullptr, -1, 3
            ));
    flooder.queue.pop();
    auto e2 = flooder.queue.top();
    ASSERT_EQ(*e2,
            pm::TentativeEvent(
                    &flooder.graph.nodes[0], 1,
                    &flooder.graph.nodes[1], 0, 5
                    )
            );
    flooder.queue.pop();
    ASSERT_EQ(flooder.graph.nodes[0].region_that_arrived->shell_area[0], &flooder.graph.nodes[0]);
    ASSERT_EQ(flooder.graph.nodes[0].distance_from_source, 0);
    ASSERT_EQ(flooder.graph.nodes[0].reached_from_source, &flooder.graph.nodes[0]);
    auto e3 = flooder.queue.top();
    ASSERT_EQ(*e3,
              pm::TentativeEvent(
                      &flooder.graph.nodes[2], 0,
                      &flooder.graph.nodes[1], 1,
                      11
                      )
              );
    flooder.queue.pop();
    auto e4 = flooder.queue.top();
    ASSERT_EQ(*e4,
              pm::TentativeEvent(
                      &flooder.graph.nodes[3], 0,
                      &flooder.graph.nodes[2], 1, 50
                      )
              );
    flooder.queue.pop();
    auto e5 = flooder.queue.top();
    auto e5_exp = pm::TentativeEvent(
            &flooder.graph.nodes[2], 1,
            &flooder.graph.nodes[3], 0,
            100
    );
    e5_exp.is_invalidated = true;
    ASSERT_EQ(*e5, e5_exp);
    flooder.queue.pop();
    auto e6 = flooder.queue.top();
    ASSERT_EQ(*e6,
              pm::TentativeEvent(
                      &flooder.graph.nodes[3], 1,
                      nullptr, -1, 1000
                      )
              );
    ASSERT_FALSE(flooder.queue.empty());
    flooder.queue.pop();
    ASSERT_TRUE(flooder.queue.empty());
}
