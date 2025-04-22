#include "pymatching/sparse_blossom/driver/detector_edge_id.h"

#include <iostream>
#include <sstream>

namespace stim_util {

DetectorEdgeId::DetectorEdgeId() : d1(0), d2(0) {
}

DetectorEdgeId::DetectorEdgeId(uint64_t d) : d1(d), d2(d) {
}

DetectorEdgeId::DetectorEdgeId(uint64_t init_d1, uint64_t init_d2) : d1(init_d1), d2(init_d2) {
    if (d2 < d1) {
        std::swap(d1, d2);
    }
}

bool DetectorEdgeId::is_boundary_edge() const {
    return d1 == d2;
}

bool DetectorEdgeId::operator==(const DetectorEdgeId &other) const {
    return d1 == other.d1 && d2 == other.d2;
}

bool DetectorEdgeId::operator!=(const DetectorEdgeId &other) const {
    return !(*this == other);
}

bool DetectorEdgeId::operator<(const DetectorEdgeId &other) const {
    return d1 != other.d1 ? d1 < other.d1 : d2 < other.d2;
}

std::ostream &operator<<(std::ostream &out, const DetectorEdgeId &edge) {
    out << "D" << edge.d1;
    if (!edge.is_boundary_edge()) {
        out << " D" << edge.d2;
    }
    return out;
}

std::string DetectorEdgeId::str() const {
    std::stringstream s;
    s << *this;
    return s.str();
}

const uint64_t *DetectorEdgeId::begin() const {
    return &d1;
}

const uint64_t *DetectorEdgeId::end() const {
    return begin() + (is_boundary_edge() ? 1 : 2);
}

/// Creates a new `DetectorEdgeId` that has detector indices offset by the given `d_offset`
/// argument.
DetectorEdgeId offset_DetectorEdgeId(const DetectorEdgeId &deid, size_t d_offset) {
    return {deid.d1 + d_offset, deid.d2 + d_offset};
}

}  // namespace stim_util
