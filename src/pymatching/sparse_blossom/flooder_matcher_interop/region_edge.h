// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PYMATCHING2_REGION_EDGE_H
#define PYMATCHING2_REGION_EDGE_H

#include <tuple>

#include "pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.h"

namespace pm {

class GraphFillRegion;

struct RegionEdge {
    GraphFillRegion* region;
    CompressedEdge edge;
    bool operator==(const RegionEdge& rhs) const;
    bool operator!=(const RegionEdge& rhs) const;
};

}  // namespace pm

#endif  // PYMATCHING2_REGION_EDGE_H
