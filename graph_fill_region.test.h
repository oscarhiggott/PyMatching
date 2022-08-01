#ifndef PYMATCHING2_GRAPH_FILL_REGION_TEST_H
#define PYMATCHING2_GRAPH_FILL_REGION_TEST_H

#include "graph_fill_region.h"
#include "graph.h"


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


#endif //PYMATCHING2_GRAPH_FILL_REGION_TEST_H
