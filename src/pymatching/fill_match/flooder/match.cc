#include "pymatching/fill_match/flooder/match.h"

bool pm::Match::operator==(const pm::Match &rhs) const {
    return region == rhs.region && edge == rhs.edge;
}

bool pm::Match::operator!=(const pm::Match &rhs) const {
    return !(rhs == *this);
}
