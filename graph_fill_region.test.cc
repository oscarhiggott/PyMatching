#include <gtest/gtest.h>
#include "graph_fill_region.h"
#include "region_edge.h"


struct GraphFillTestData {
    std::vector<pm::DetectorNode> detectors;
    explicit GraphFillTestData(int num_elements) {
        detectors.resize(num_elements);
    };
    pm::RegionEdge b(int loc_from, int loc_to, std::vector<pm::RegionEdge> edges,
                     bool root = false) {
        auto r = new pm::GraphFillRegion();
        r->blossom_children = std::move(edges);
        for (auto c : r->blossom_children)
            c.region->blossom_parent = r;
        auto ce = root
                ?  pm::CompressedEdge()
                : pm::CompressedEdge(&detectors[loc_from], &detectors[loc_to], 0);
        auto region_edge = pm::RegionEdge(
                r,
                ce);
        return region_edge;
    }
};


TEST(GraphFillRegion, BlossomBuilder) {
    GraphFillTestData x(4);

    auto t1 = x.b(-1, -1, {
        x.b(1, 2, {}, false), x.b(2, 3, {}, false), x.b(3, 1, {}, false)
        }, true);
    auto t2 = x.b(-1, -1, {
            x.b(1, 2, {}, false), x.b(2, 3, {}, false), x.b(3, 1, {}, false)
    }, true);
    ASSERT_EQ(*t1.region, *t2.region);
}


TEST(GraphFillRegion, TopRegion) {
    GraphFillTestData x(10);
    auto t = x.b(-1, -1, {
        x.b(0, 1, {}, false),
        x.b(1, 2, {}, false),
        x.b(2, 0, {
            x.b(3, 4, {}, false),
            x.b(4, 5, {}, false),
            x.b(5, 3, {}, false)
            }, false)
        }, true);

    ASSERT_EQ(t.region, t.region->blossom_children[2].region->blossom_children[0].region->top_region());
}
