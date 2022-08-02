#ifndef PYMATCHING2_COMPRESSED_EDGE_H
#define PYMATCHING2_COMPRESSED_EDGE_H

#include "graph.h"

namespace pm {


    struct CompressedEdge{
        DetectorNode* loc_from;
        DetectorNode* loc_to;
        obs_int obs_mask;

        CompressedEdge reversed() const;

        bool operator==(const CompressedEdge &rhs) const;

        bool operator!=(const CompressedEdge &rhs) const;

        CompressedEdge merged_with(const CompressedEdge& other) const;
        CompressedEdge() = default;
        CompressedEdge(DetectorNode* loc_from, DetectorNode* loc_to, obs_int obs_mask);
    };

    inline CompressedEdge::CompressedEdge(DetectorNode* loc_from, DetectorNode* loc_to, obs_int obs_mask)
        : loc_from(loc_from), loc_to(loc_to), obs_mask(obs_mask) {}

    inline CompressedEdge CompressedEdge::reversed() const {
        return {loc_to, loc_from, obs_mask};
    }

    inline CompressedEdge CompressedEdge::merged_with(const CompressedEdge &other) const {
        return {loc_from, other.loc_to, obs_mask ^ other.obs_mask};
    }

}

#endif //PYMATCHING2_COMPRESSED_EDGE_H
