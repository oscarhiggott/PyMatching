#ifndef PYMATCHING2_GRAPH_FILL_REGION_TEST_H
#define PYMATCHING2_GRAPH_FILL_REGION_TEST_H

#include "pymatching/fill_match/flooder/graph.h"
#include "pymatching/fill_match/flooder/graph_fill_region.h"
#include "pymatching/fill_match/arena.h"

namespace pm {
struct GraphFillTestData {
    std::vector<DetectorNode> detectors;
    Arena<GraphFillRegion> arena;
    explicit GraphFillTestData(int num_elements) {
        detectors.resize(num_elements);
    };
    RegionEdge b(int loc_from, int loc_to, std::vector<RegionEdge> edges, bool root = false) {
        auto r = arena.alloc_default_constructed();
        r->blossom_children = std::move(edges);
        for (auto c : r->blossom_children) {
            c.region->wrap_into_blossom(r);
        }
        auto ce = root ? CompressedEdge{nullptr, nullptr, 0} : CompressedEdge{&detectors[loc_from], &detectors[loc_to], 0};
        auto region_edge = RegionEdge{r, ce};
        return region_edge;
    }
};

}  // namespace pm

#endif  // PYMATCHING2_GRAPH_FILL_REGION_TEST_H
