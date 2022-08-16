#include "pymatching/fill_match/flooder_matcher_interop/compressed_edge.h"

bool pm::CompressedEdge::operator==(const pm::CompressedEdge &rhs) const {
    return loc_from == rhs.loc_from && loc_to == rhs.loc_to && obs_mask == rhs.obs_mask;
}

bool pm::CompressedEdge::operator!=(const pm::CompressedEdge &rhs) const {
    return !(rhs == *this);
}

std::string pm::CompressedEdge::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

std::ostream &pm::operator<<(std::ostream &out, const pm::CompressedEdge &edge) {
    out << "CompressedEdge{.obs_mask=" << edge.obs_mask;
    out << ", .loc_from=" << edge.loc_from;
    out << ", .loc_to=" << edge.loc_to;
    out << "}";
    return out;
}
