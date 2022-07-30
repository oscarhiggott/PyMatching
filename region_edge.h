#ifndef PYMATCHING2_REGION_EDGE_H
#define PYMATCHING2_REGION_EDGE_H

#include "compressed_edge.h"
#include <tuple>

namespace pm {

    class GraphFillRegion;

// Instead blossom children just has vector of regions and edges. Write separate functions that take both

    struct RegionEdge{
        GraphFillRegion* region;
        CompressedEdge edge;
        RegionEdge();
        RegionEdge(GraphFillRegion* region, CompressedEdge edge);

        bool operator==(const RegionEdge &rhs) const;

        bool operator!=(const RegionEdge &rhs) const;
    };


//    class RegionPath {
//    public:
//        std::vector<RegionEdge> edges;
//        std::tuple<RegionPath,RegionPath> split_between_regions(GraphFillRegion* start_region, GraphFillRegion* end_region);
//        RegionPath split_at_region(GraphFillRegion* region);
//        std::vector<GraphFillRegion*> pairs_matched();
//    };


}

#endif //PYMATCHING2_REGION_EDGE_H
