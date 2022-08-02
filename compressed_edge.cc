#include "compressed_edge.h"


bool pm::CompressedEdge::operator==(const pm::CompressedEdge &rhs) const {
    return loc_from == rhs.loc_from &&
           loc_to == rhs.loc_to &&
           obs_mask == rhs.obs_mask;
}

bool pm::CompressedEdge::operator!=(const pm::CompressedEdge &rhs) const {
    return !(rhs == *this);
}
