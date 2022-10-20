#include "pymatching/sparse_blossom/flooder_matcher_interop/region_edge.h"

bool pm::RegionEdge::operator==(const pm::RegionEdge &rhs) const {
    return region == rhs.region && edge == rhs.edge;
}

bool pm::RegionEdge::operator!=(const pm::RegionEdge &rhs) const {
    return !(rhs == *this);
}
