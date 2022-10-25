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

#ifndef PYMATCHING2_SEARCH_GRAPH_H
#define PYMATCHING2_SEARCH_GRAPH_H

#include "pymatching/sparse_blossom/search/search_detector_node.h"
#include "pymatching/sparse_blossom/tracker/queued_event_tracker.h"

namespace pm {

/// The edge on which two search regions collided
struct SearchGraphEdge {
    SearchDetectorNode* detector_node;
    size_t neighbor_index;
};

class SearchGraph {
   public:
    std::vector<SearchDetectorNode> nodes;
    size_t num_nodes;

    SearchGraph();
    explicit SearchGraph(size_t num_nodes);
    SearchGraph(SearchGraph&& graph) noexcept;
    void add_edge(size_t u, size_t v, signed_weight_int weight, const std::vector<size_t>& observables);
    void add_boundary_edge(size_t u, signed_weight_int weight, const std::vector<size_t>& observables);
};
}  // namespace pm

#endif  // PYMATCHING2_SEARCH_GRAPH_H
