#include "pymatching/sparse_blossom/flooder/match.h"

bool pm::Match::operator==(const pm::Match &rhs) const {
    return region == rhs.region && edge == rhs.edge;
}

bool pm::Match::operator!=(const pm::Match &rhs) const {
    return !(rhs == *this);
}
