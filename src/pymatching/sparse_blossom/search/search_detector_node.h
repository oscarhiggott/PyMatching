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

#ifndef PYMATCHING2_SEARCH_DETECTOR_NODE_H
#define PYMATCHING2_SEARCH_DETECTOR_NODE_H

#include "pymatching/sparse_blossom/tracker/queued_event_tracker.h"

namespace pm {

class SearchDetectorNode {
   public:
    SearchDetectorNode() : reached_from_source(nullptr), index_of_predecessor(SIZE_MAX), distance_from_source(0) {
    }

    /// The SearchDetectorNode that this node was reached from in the Dijkstra search
    SearchDetectorNode *reached_from_source;

    /// `index_of_predecessor` is the index in `neighbors` of the neighboring detector node that this node was
    /// reached from in the Dijkstra search.
    size_t index_of_predecessor;

    /// The distance from the search source detector node that this node was reached from.
    pm::cumulative_time_int distance_from_source;

    /// Manages the next "look at me" event for the node
    QueuedEventTracker node_event_tracker;

    /// == Permanent fields used to define the structure of the graph. ==
    std::vector<SearchDetectorNode *> neighbors;  /// The node's neighbors.
    std::vector<weight_int> neighbor_weights;     /// Distance crossed by the edge to each neighbor.
    std::vector<std::vector<size_t>>
        neighbor_observable_indices;  /// Indices of observables crossed by the edge to each neighbor.

    size_t index_of_neighbor(SearchDetectorNode *target) const;

    void reset();
};

}  // namespace pm

#endif  // PYMATCHING2_SEARCH_DETECTOR_NODE_H
