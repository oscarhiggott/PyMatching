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

#ifndef PYMATCHING2_GRAPH_FILL_REGION_TEST_H
#define PYMATCHING2_GRAPH_FILL_REGION_TEST_H

#include "pymatching/sparse_blossom/arena.h"
#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"

namespace pm {
struct GraphFillTestData {
    std::vector<DetectorNode> detectors;
    Arena<GraphFillRegion> arena;
    explicit GraphFillTestData(int num_elements) {
        detectors.resize(num_elements);
    };
    RegionEdge b(int loc_from, int loc_to, std::vector<RegionEdge> edges, bool root = false) {
        auto r = arena.alloc_default_constructed();
        r->blossom_children = std::move(edges);
        for (auto c : r->blossom_children) {
            c.region->wrap_into_blossom(r);
        }
        auto ce =
            root ? CompressedEdge{nullptr, nullptr, 0} : CompressedEdge{&detectors[loc_from], &detectors[loc_to], 0};
        auto region_edge = RegionEdge{r, ce};
        return region_edge;
    }
};

}  // namespace pm

#endif  // PYMATCHING2_GRAPH_FILL_REGION_TEST_H
