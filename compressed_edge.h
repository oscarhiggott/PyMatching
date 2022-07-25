#ifndef PYMATCHING2_COMPRESSED_EDGE_H
#define PYMATCHING2_COMPRESSED_EDGE_H

#include "graph.h"

namespace pm {


    struct CompressedEdge{
        DetectorNode* loc_from;
        DetectorNode* loc_to;
        obs_int obs_mask;

        CompressedEdge reversed();
        CompressedEdge merged_with(CompressedEdge other);
    };

}

#endif //PYMATCHING2_COMPRESSED_EDGE_H
