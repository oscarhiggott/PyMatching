#include "pymatching/fill_match/flooder/graph_fill_region.h"

#include <gtest/gtest.h>

#include "pymatching/fill_match/flooder/graph_fill_region.test.h"
#include "pymatching/fill_match/flooder_matcher_interop/region_edge.h"

using namespace pm;

TEST(GraphFillRegion, BlossomBuilder) {
    GraphFillTestData x(4);

    auto t1 = x.b(-1, -1, {x.b(1, 2, {}, false), x.b(2, 3, {}, false), x.b(3, 1, {}, false)}, true);
    auto t2 = x.b(-1, -1, {x.b(1, 2, {}, false), x.b(2, 3, {}, false), x.b(3, 1, {}, false)}, true);
    ASSERT_EQ(*t1.region, *t2.region);
}

TEST(GraphFillRegion, TopRegion) {
    GraphFillTestData x(10);
    auto t =
        x.b(-1,
            -1,
            {x.b(0, 1, {}, false),
             x.b(1, 2, {}, false),
             x.b(2, 0, {x.b(3, 4, {}, false), x.b(4, 5, {}, false), x.b(5, 3, {}, false)}, false)},
            true);

    ASSERT_EQ(t.region, t.region->blossom_children[2].region->blossom_children[0].region->blossom_parent_top);
    ASSERT_EQ(t.region, t.region->blossom_children[1].region->blossom_parent_top);
}

TEST(GraphFillRegion, AddMatch) {
    GraphFillRegion r1;
    GraphFillRegion r2;
    std::vector<DetectorNode> ds(2);
    r1.add_match(&r2, CompressedEdge{&ds[0], &ds[1], 5});
    ASSERT_EQ(r1.match.region, &r2);
    ASSERT_EQ(r1.match.edge, r2.match.edge.reversed());
    ASSERT_EQ(r1.match.edge, (CompressedEdge{&ds[0], &ds[1], 5}));
}
