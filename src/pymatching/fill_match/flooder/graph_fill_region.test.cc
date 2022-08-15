#include "pymatching/fill_match/flooder/graph_fill_region.h"

#include <gtest/gtest.h>

#include "pymatching/fill_match/flooder/graph_fill_region.test.h"
#include "pymatching/fill_match/matcher/region_edge.h"

void pm::delete_graph_fill_region(pm::GraphFillRegion* region) {
    for (auto& child : region->blossom_children)
        delete_graph_fill_region(child.region);
    delete region;
}

TEST(GraphFillRegion, BlossomBuilder) {
    pm::GraphFillTestData x(4);

    auto t1 = x.b(-1, -1, {x.b(1, 2, {}, false), x.b(2, 3, {}, false), x.b(3, 1, {}, false)}, true);
    auto t2 = x.b(-1, -1, {x.b(1, 2, {}, false), x.b(2, 3, {}, false), x.b(3, 1, {}, false)}, true);
    ASSERT_EQ(*t1.region, *t2.region);
    pm::delete_graph_fill_region(t1.region);
    pm::delete_graph_fill_region(t2.region);
}

TEST(GraphFillRegion, TopRegion) {
    pm::GraphFillTestData x(10);
    auto t =
        x.b(-1,
            -1,
            {x.b(0, 1, {}, false),
             x.b(1, 2, {}, false),
             x.b(2, 0, {x.b(3, 4, {}, false), x.b(4, 5, {}, false), x.b(5, 3, {}, false)}, false)},
            true);

    ASSERT_EQ(t.region, t.region->blossom_children[2].region->blossom_children[0].region->top_region());
    ASSERT_EQ(t.region, t.region->blossom_children[1].region->top_region());
    pm::delete_graph_fill_region(t.region);
}

TEST(GraphFillRegion, AddMatch) {
    pm::GraphFillRegion r1;
    pm::GraphFillRegion r2;
    std::vector<pm::DetectorNode> ds(2);
    r1.add_match(&r2, pm::CompressedEdge(&ds[0], &ds[1], 5));
    ASSERT_EQ(r1.match.region, &r2);
    ASSERT_EQ(r1.match.edge, r2.match.edge.reversed());
    ASSERT_EQ(r1.match.edge, pm::CompressedEdge(&ds[0], &ds[1], 5));
}
