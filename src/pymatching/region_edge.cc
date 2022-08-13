#include "pymatching/region_edge.h"

pm::RegionEdge::RegionEdge(pm::GraphFillRegion *region, CompressedEdge edge) : region(region), edge(edge) {
}

pm::RegionEdge::RegionEdge() : region(nullptr), edge(CompressedEdge()) {
}

bool pm::RegionEdge::operator==(const pm::RegionEdge &rhs) const {
    return region == rhs.region && edge == rhs.edge;
}

bool pm::RegionEdge::operator!=(const pm::RegionEdge &rhs) const {
    return !(rhs == *this);
}
