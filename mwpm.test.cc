#include <gtest/gtest.h>
#include "mwpm.h"


struct MwpmAltTreeTestData {
    std::vector<pm::DetectorNode> nodes;
    explicit MwpmAltTreeTestData(std::vector<pm::DetectorNode>& nodes);
    pm::AltTreeEdge t(
            std::vector<pm::AltTreeEdge> children,
            size_t inner_region_id,
            size_t outer_region_id,
            bool root = false
    );
};

MwpmAltTreeTestData::MwpmAltTreeTestData(std::vector<pm::DetectorNode>& nodes) : nodes(std::move(nodes)) {}

pm::AltTreeEdge
MwpmAltTreeTestData::t(std::vector<pm::AltTreeEdge> children, size_t inner_region_id, size_t outer_region_id, bool root) {
    pm::AltTreeNode* node;
    pm::CompressedEdge parent_ce(nullptr, nullptr, 0);
    if (root) {
        node = new pm::AltTreeNode(nodes[outer_region_id].region_that_arrived);
    } else {
        node = new pm::AltTreeNode(
                nodes[inner_region_id].region_that_arrived,
                nodes[outer_region_id].region_that_arrived,
                {
                        &nodes[inner_region_id],
                        &nodes[outer_region_id],
                        0
                }
        );
    }
    auto edge = pm::AltTreeEdge(node, parent_ce);
    for (auto& child : children) {
        child.edge.loc_from = &nodes[outer_region_id];
        child.edge.loc_to = child.alt_tree_node->inner_to_outer_edge.loc_from;
        edge.alt_tree_node->add_child(child);
    }
    return edge;
}

TEST(Mwpm, BlossomCreatedThenShattered) {
    auto g = pm::Graph(10);
    g.add_edge(0, 1, 10, 1);
    g.add_edge(1, 4, 20, 2);
    g.add_edge(4, 3, 20, 3);
    g.add_edge(3, 2, 12, 4);
    g.add_edge(0, 2, 16, 5);
    g.add_edge(4, 5, 50, 6);
    g.add_edge(2, 6, 100, 7);
    g.add_boundary_edge(5, 36, 8);
    auto flooder = pm::GraphFlooder(g);
    auto mwpm = pm::Mwpm(flooder);
    for (size_t i = 0; i < 7; i++)
        mwpm.flooder.create_region(&mwpm.flooder.graph.nodes[i]);
    auto& ns = mwpm.flooder.graph.nodes;
    auto e1 = mwpm.flooder.next_event();
    mwpm.process_event(e1);
    auto e2 = mwpm.flooder.next_event();
    mwpm.process_event(e2);
    auto e3 = mwpm.flooder.next_event();
    mwpm.process_event(e3);
    auto e4 = mwpm.flooder.next_event();
    mwpm.process_event(e4);
    auto e5 = mwpm.flooder.next_event();
    ASSERT_EQ(e5, pm::MwpmEvent(
            mwpm.flooder.graph.nodes[0].region_that_arrived,
            mwpm.flooder.graph.nodes[2].region_that_arrived,
            {&mwpm.flooder.graph.nodes[0], &mwpm.flooder.graph.nodes[2], 5}
            ));
    mwpm.process_event(e5);
    // Region 5 matches with blossom containing {0, 1, 4, 3, 2}
    auto e6 = mwpm.flooder.next_event();
    mwpm.process_event(e6);
    auto blossom = mwpm.flooder.graph.nodes[5].region_that_arrived->match.region;
    ASSERT_EQ(blossom->blossom_children.size(), 5);
    std::vector<pm::RegionEdge> expected_blossom_children = {
            {ns[2].region_that_arrived, {&ns[2], &ns[3], 4}},
            {ns[3].region_that_arrived, {&ns[3], &ns[4], 3}},
            {ns[4].region_that_arrived, {&ns[4], &ns[1], 2}},
            {ns[1].region_that_arrived, {&ns[1], &ns[0], 1}},
            {ns[0].region_that_arrived, {&ns[0], &ns[2], 5}},
    };
    ASSERT_EQ(blossom->blossom_children, expected_blossom_children);
    ASSERT_EQ(blossom->radius, pm::Varying(8 << 2));
    ASSERT_EQ(mwpm.flooder.time, 25);
    // Region 6 collides with matched blossom
    auto e7 = mwpm.flooder.next_event();
    mwpm.process_event(e7);
    ASSERT_EQ(mwpm.flooder.time, 83);
    ASSERT_EQ(
            mwpm.flooder.graph.nodes[6].region_that_arrived->alt_tree_node->children[0].alt_tree_node->inner_region,
            blossom
              );
    // Blossom shatters
    auto e8 = mwpm.flooder.next_event();
    ASSERT_EQ(mwpm.flooder.time, 91);
    ASSERT_EQ(e8.event_type, pm::BLOSSOM_SHATTER);
    mwpm.process_event(e8);
    ASSERT_EQ(ns[0].region_that_arrived->match.region,
              ns[1].region_that_arrived);
    ASSERT_EQ(ns[1].region_that_arrived->match.region,
              ns[0].region_that_arrived);
    auto e9 = mwpm.flooder.next_event();
    ASSERT_EQ(mwpm.flooder.time, 94);
    ASSERT_EQ(e9.event_type, pm::REGION_HIT_BOUNDARY);
    mwpm.process_event(e9);
    auto e10 = mwpm.flooder.next_event();
    ASSERT_EQ(e10.event_type, pm::NO_EVENT);
    ASSERT_EQ(ns[3].region_that_arrived->blossom_parent, nullptr);
    auto regions_matched = [&ns](size_t i, size_t j) {
        bool regions_correct = ns[i].region_that_arrived->match.region == ns[j].region_that_arrived &&
                                ns[j].region_that_arrived->match.region == ns[i].region_that_arrived;
        bool edge_locs_correct = ns[i].region_that_arrived->match.edge.loc_from == &ns[i] &&
                                  ns[i].region_that_arrived->match.edge.loc_to == &ns[j] &&
                                  ns[j].region_that_arrived->match.edge.loc_from == &ns[j] &&
                                  ns[j].region_that_arrived->match.edge.loc_to == &ns[i];
        return regions_correct && edge_locs_correct;
    };
    ASSERT_TRUE(regions_matched(4, 3));
    ASSERT_TRUE(regions_matched(2, 6));
    ASSERT_TRUE(regions_matched(1, 0));
    ASSERT_EQ(ns[5].region_that_arrived->match, pm::Match(
            nullptr, {&ns[5], nullptr, 8}
            ));
}

TEST(Mwpm, BlossomShatterDrivenWithoutFlooder) {
    size_t n = 10;
    pm::Graph g(n+3);
    auto flooder = pm::GraphFlooder(g);
    auto mwpm = pm::Mwpm(flooder);
    auto& ns = mwpm.flooder.graph.nodes;
    for (size_t i = 0; i < n + 3; i++)
        mwpm.flooder.create_region(&ns[i]);

    auto rhr = [&ns](size_t i, size_t j) {
        return pm::MwpmEvent(
                ns[i].region_that_arrived,
                ns[j].region_that_arrived,
                {&ns[i], &ns[j], 0}
                );
    };
    // Pair up
    for (size_t i = 0; i < n; i += 2)
        mwpm.process_event(rhr(i, i+1));
    // Form an alternating path
    for (size_t i = 1; i < n; i += 2)
        mwpm.process_event(rhr(n-i, n-i+1));
    // Close the path into a blossom
    mwpm.process_event(rhr(0, n));
    auto blossom_id = n + 3;
    auto blossom = ns[0].region_that_arrived->blossom_parent;
    // Make the blossom become an inner node
    mwpm.process_event(
            {
                ns[n+1].region_that_arrived,
                blossom,
                {&ns[n+1], &ns[5], 0}
            }
            );
    mwpm.process_event(
            {
                    ns[n+2].region_that_arrived,
                    blossom,
                    {&ns[n+2], &ns[10], 0}
            }
    );
    ASSERT_EQ(blossom->blossom_children.size(), 11);
    mwpm.process_event(
            {
              blossom,
             ns[10].region_that_arrived,
             ns[5].region_that_arrived
            }
            );
    auto regions_matched = [&ns](size_t i, size_t j) {
        bool regions_correct = ns[i].region_that_arrived->match.region == ns[j].region_that_arrived &&
                               ns[j].region_that_arrived->match.region == ns[i].region_that_arrived;
        bool edge_locs_correct = ns[i].region_that_arrived->match.edge.loc_from == &ns[i] &&
                                 ns[i].region_that_arrived->match.edge.loc_to == &ns[j] &&
                                 ns[j].region_that_arrived->match.edge.loc_from == &ns[j] &&
                                 ns[j].region_that_arrived->match.edge.loc_to == &ns[i];
        return regions_correct && edge_locs_correct;
    };
    ASSERT_TRUE(regions_matched(8, 9));
    ASSERT_TRUE(regions_matched(6, 7));
    auto actual_tree = ns[12].region_that_arrived->alt_tree_node;
    ASSERT_EQ(actual_tree->parent.alt_tree_node, nullptr);

    MwpmAltTreeTestData d(ns);
    auto expected_tree = d.t(
            {
                    d.t({
                        d.t({
                            d.t({
                                d.t({}, 5, 11)
                                }, 3, 4)
                            }, 1, 2)
                        }, 10, 0)
                },
            -1, 12, true
            ).alt_tree_node;
    ASSERT_EQ(*actual_tree, *expected_tree);
}
