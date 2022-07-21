#ifndef PYMATCHING2_REGION_PATH_H
#define PYMATCHING2_REGION_PATH_H

#include "compressed_edge.h"
#include <tuple>


class GraphFillRegion;

// Instead blossom children just has vector of regions and edges. Write separate functions that take both

struct RegionEdge{
    GraphFillRegion* region;
    CompressedEdge edge;
};


class RegionPath {
public:
    std::vector<RegionEdge> edges;
    std::tuple<RegionPath,RegionPath> split_between_regions(GraphFillRegion* start_region, GraphFillRegion* end_region);
    RegionPath split_at_region(GraphFillRegion* region);
    std::vector<GraphFillRegion*> pairs_matched();
};


#endif //PYMATCHING2_REGION_PATH_H
