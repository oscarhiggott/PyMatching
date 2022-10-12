#ifndef PYMATCHING2_COMPRESSED_EDGE_H
#define PYMATCHING2_COMPRESSED_EDGE_H

#include "pymatching/fill_match/ints.h"

namespace pm {

struct DetectorNode;

/// A compressed edge is a summary of a path between two detection events.
/// Specifically, it tracks which observables the path has crossed.
///
/// The matching algorithm tracks compressed edges, instead of entire
/// paths, because this is more space efficient and more time efficient
/// and because the full path information is not actually needed in order
/// to predict whether or not observables were flipped (which is the
/// minimalist goal of the decoder). (When full path information is
/// requested, it can be reconstructed from the compressed edges
/// after matching has finished by using Dijkstra's algorithm.)
struct CompressedEdge {
    /// The detection event where the path starts.
    DetectorNode* loc_from;
    /// The detection event where the path ends.
    DetectorNode* loc_to;
    /// A bit mask of which observables were crossed an odd number of
    /// times by the path. The bit 1<<K is set if the K'th observable
    /// was crossed an odd number of times.
    obs_int obs_mask;

    CompressedEdge reversed() const;

    bool operator==(const CompressedEdge& rhs) const;
    bool operator!=(const CompressedEdge& rhs) const;

    CompressedEdge merged_with(const CompressedEdge& other) const;

    std::string str() const;

    inline void clear() {
        loc_from = nullptr;
        loc_to = nullptr;
        obs_mask = 0;
    }
};

inline CompressedEdge CompressedEdge::reversed() const {
    return {loc_to, loc_from, obs_mask};
}

inline CompressedEdge CompressedEdge::merged_with(const CompressedEdge& other) const {
    return {loc_from, other.loc_to, obs_mask ^ other.obs_mask};
}

std::ostream& operator<<(std::ostream& out, const pm::CompressedEdge& edge);

}  // namespace pm

#endif  // PYMATCHING2_COMPRESSED_EDGE_H
