#ifndef PYMATCHING2_REGION_EDGE_H
#define PYMATCHING2_REGION_EDGE_H

#include "compressed_edge.h"
#include <tuple>

namespace pm {

    class GraphFillRegion;

    struct RegionEdge{
        GraphFillRegion* region;
        CompressedEdge edge;
        RegionEdge();
        RegionEdge(GraphFillRegion* region, CompressedEdge edge);

        bool operator==(const RegionEdge &rhs) const;

        bool operator!=(const RegionEdge &rhs) const;
    };

}

#endif //PYMATCHING2_REGION_EDGE_H
