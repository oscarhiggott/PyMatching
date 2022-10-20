#ifndef PYMATCHING_FILL_MATCH_MATCH_H
#define PYMATCHING_FILL_MATCH_MATCH_H

#include <vector>

#include "pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.h"

namespace pm {

class GraphFillRegion;

/// A Match is a partnered graph fill region.
///
/// When one graph fill region is matched to another, it's guaranteed that (if the matching
/// process terminated) two of their detection events will be matched to each other. The
/// specific detection events to match are identified by `edge.loc_from` and `edge.loc_to`.
struct Match {
    /// The region being matched to (or nullptr if matched to the boundary).
    pm::GraphFillRegion* region;
    /// A summary of the low-level path from this region to that region, including the
    /// start/end detection events and the observables that were crossed.
    pm::CompressedEdge edge;
    bool operator==(const Match& rhs) const;
    bool operator!=(const Match& rhs) const;

    inline void clear() {
        region = nullptr;
        edge.clear();
    }
};

}  // namespace pm

#endif  // PYMATCHING_FILL_MATCH_MATCH_H
