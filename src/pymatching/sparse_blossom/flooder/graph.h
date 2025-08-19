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

#ifndef PYMATCHING2_GRAPH_H
#define PYMATCHING2_GRAPH_H

#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#include "pymatching/sparse_blossom/driver/implied_weights.h"
#include "pymatching/sparse_blossom/flooder/detector_node.h"
#include "pymatching/sparse_blossom/tracker/flood_check_event.h"

namespace pm {

struct GraphFillRegion;

struct PreviousWeight {
    weight_int* ptr;
    weight_int val;
    PreviousWeight(weight_int* ptr, weight_int val) : ptr(ptr), val(val) {
    }
};

/// A collection of detector nodes. It's expected that all detector nodes in the graph
/// will only refer to other detector nodes within the same graph.
class MatchingGraph {
   public:
    std::vector<DetectorNode> nodes;
    /// These are the detection events that would occur if an error occurred on every edge with a negative weight
    std::set<size_t> negative_weight_detection_events_set;
    /// These are the observables that would be flipped if an error occurred on every edge with a negative weight
    std::set<size_t> negative_weight_observables_set;
    /// The sum of the negative edge weights. This number is negative, rather than the absolute value.
    pm::total_weight_int negative_weight_sum;
    /// is_user_graph_boundary_node is only filled if MatchingGraph is constructed from a UserGraph
    /// is_user_graph_boundary_node[i] is true if i is a boudary node in the UserGraph
    /// This vector is used to check that a detection event is not added to a node marked as a boundary node
    /// in the UserGraph (as such an event cannot be matched, and will raise an error).
    std::vector<bool> is_user_graph_boundary_node;
    size_t num_nodes;
    size_t num_observables;
    /// This is the normalising constant that the edge weights were multiplied by when converting from floats to
    /// 16-bit ints.
    double normalising_constant;
    std::vector<PreviousWeight> previous_weights;
    std::map<size_t, std::vector<std::vector<pm::ImpliedWeightUnconverted>>> edges_to_implied_weights_unconverted;
    // Tracks whether the MatchingGraph was loaded from a DEM without enable_correlations. This is used to
    // alert a user if they try to decode with enable_correlations=true, but forgot to load from the
    // dem with enable_correlations=true.
    bool loaded_from_dem_without_correlations = false;

    MatchingGraph();
    MatchingGraph(size_t num_nodes, size_t num_observables);
    MatchingGraph(size_t num_nodes, size_t num_observables, double normalising_constant);
    MatchingGraph(MatchingGraph&& graph) noexcept;
    void add_edge(
        size_t u,
        size_t v,
        signed_weight_int weight,
        const std::vector<size_t>& observables,
        const std::vector<pm::ImpliedWeightUnconverted>& implied_weights_for_other_edges = {});
    void add_boundary_edge(
        size_t u,
        signed_weight_int weight,
        const std::vector<size_t>& observables,
        const std::vector<pm::ImpliedWeightUnconverted>& implied_weights_for_other_edges = {});
    void update_negative_weight_observables(const std::vector<size_t>& observables);
    void update_negative_weight_detection_events(size_t node_id);
    void convert_implied_weights(double normalising_constant);

    void undo_reweights();
    void reweight(std::vector<ImpliedWeight>& implied_weights);
    void reweight_for_edge(const int64_t& u, const int64_t& v);
    void reweight_for_edges(const std::vector<int64_t>& edges);
};

void apply_reweights(
    std::vector<std::tuple<pm::weight_int*, pm::weight_int*, pm::weight_int>>& implied_weights,
    std::vector<PreviousWeight>& previous_weights);

inline void apply_reweights(
    std::vector<ImpliedWeight>& implied_weights, std::vector<PreviousWeight>& previous_weights) {
    for (auto& [edge0_ptr, edge1_ptr, new_weight] : implied_weights) {
        // Only reweight if the new weight is lower than the current weight
        if (new_weight < *edge0_ptr) {
            previous_weights.emplace_back(edge0_ptr, *edge0_ptr);
            *edge0_ptr = new_weight;
            if (edge1_ptr != nullptr) {
                // We already know new_weight < *edge1_ptr, since *edge0_ptr ==
                // *edge1_ptr
                previous_weights.emplace_back(edge1_ptr, *edge1_ptr);
                *edge1_ptr = new_weight;
            }
        }
    }
}

inline void MatchingGraph::reweight(std::vector<ImpliedWeight>& implied_weights) {
    apply_reweights(implied_weights, previous_weights);
}

}  // namespace pm

#endif  // PYMATCHING2_GRAPH_H
