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

#ifndef PYMATCHING2_SEARCH_FLOODER_H
#define PYMATCHING2_SEARCH_FLOODER_H

#include "pymatching/sparse_blossom/search/search_graph.h"
#include "pymatching/sparse_blossom/tracker/radix_heap_queue.h"

namespace pm {

enum TargetType : uint8_t { DETECTOR_NODE, BOUNDARY, NO_TARGET };

class SearchFlooder {
   public:
    SearchFlooder();
    explicit SearchFlooder(SearchGraph graph);
    SearchFlooder(SearchFlooder&& other) noexcept;
    SearchGraph graph;
    pm::radix_heap_queue<false> queue;
    /// The reached_nodes are the nodes that need to be reset after each search completes
    std::vector<SearchDetectorNode*> reached_nodes;
    /// The type of target for the search from a detection event, either another detection event or the boundary.
    TargetType target_type;
    void reschedule_events_at_search_detector_node(SearchDetectorNode& detector_node);
    std::pair<size_t, cumulative_time_int> find_next_event_at_node_returning_neighbor_index_and_time(
        const SearchDetectorNode& detector_node) const;
    void do_search_starting_at_empty_search_detector_node(SearchDetectorNode* src);
    void do_search_exploring_empty_detector_node(SearchDetectorNode& empty_node, size_t empty_to_from_index);
    SearchGraphEdge do_look_at_node_event(SearchDetectorNode& node);
    SearchGraphEdge run_until_collision(SearchDetectorNode* src, SearchDetectorNode* dst);
    template <typename Callable>
    void iter_edges_on_path_traced_back_from_node(SearchDetectorNode* detector_node, Callable handle_edge);
    template <typename Callable>
    void iter_edges_tracing_back_from_collision_edge(const SearchGraphEdge& collision_edge, Callable handle_edge);
    // Visits the edges on the shortest path between src to dst, calling handle_edge on each SearchGraphEdge.
    // This method is not guaranteed to visit nodes in order from src to dst. It traces back paths from
    // the collision point of the two search regions to src and dst.
    template <typename Callable>
    void iter_edges_on_shortest_path_from_middle(size_t src, size_t dst, Callable handle_edge);
    template <typename Callable>
    void reverse_path_and_handle_edges(const std::vector<SearchGraphEdge>& edges, Callable handle_edge);
    // Visits the edges on the shortest path from src to dst, calling handle_edge on each SearchGraphEdge.
    // This method is slightly slower than iter_edges_on_shortest_path_from_middle (since it uses a stack), but
    // guarantees that edges are visited in order from src to dst.
    template <typename Callable>
    void iter_edges_on_shortest_path_from_source(size_t src, size_t dst, Callable handle_edge);
    void reset_graph();
    void reset();
};

template <typename Callable>
void SearchFlooder::iter_edges_on_path_traced_back_from_node(SearchDetectorNode* detector_node, Callable handle_edge) {
    auto current_node = detector_node;
    while (current_node->index_of_predecessor != SIZE_MAX) {
        auto pred_idx = current_node->index_of_predecessor;
        SearchGraphEdge edge = {current_node, pred_idx};
        handle_edge(edge);
        current_node = current_node->neighbors[pred_idx];
    }
}

template <typename Callable>
void SearchFlooder::iter_edges_tracing_back_from_collision_edge(
    const SearchGraphEdge& collision_edge, Callable handle_edge) {
    iter_edges_on_path_traced_back_from_node(collision_edge.detector_node, handle_edge);
    auto other_node = collision_edge.detector_node->neighbors[collision_edge.neighbor_index];
    if (other_node)
        iter_edges_on_path_traced_back_from_node(other_node, handle_edge);
    handle_edge(collision_edge);
}

template <typename Callable>
void SearchFlooder::iter_edges_on_shortest_path_from_middle(size_t src, size_t dst, Callable handle_edge) {
    SearchDetectorNode* loc_to_ptr = dst == SIZE_MAX ? nullptr : &graph.nodes[dst];
    auto collision_edge = run_until_collision(&graph.nodes[src], loc_to_ptr);
    iter_edges_tracing_back_from_collision_edge(collision_edge, handle_edge);
    reset();
}

template <typename Callable>
void SearchFlooder::reverse_path_and_handle_edges(const std::vector<SearchGraphEdge>& edges, Callable handle_edge) {
    size_t n = edges.size();
    for (size_t i = 0; i < n; i++) {
        auto e = edges[n - 1 - i];
        SearchDetectorNode* n_from = e.detector_node->neighbors[e.neighbor_index];
        SearchGraphEdge e_rev = {n_from, n_from->index_of_neighbor(e.detector_node)};
        handle_edge(e_rev);
    }
}

template <typename Callable>
void SearchFlooder::iter_edges_on_shortest_path_from_source(size_t src, size_t dst, Callable handle_edge) {
    SearchDetectorNode* loc_from_ptr = &graph.nodes[src];
    SearchDetectorNode* loc_to_ptr = dst == SIZE_MAX ? nullptr : &graph.nodes[dst];
    auto collision_edge = run_until_collision(loc_from_ptr, loc_to_ptr);

    std::vector<SearchGraphEdge> path_edges_1;
    iter_edges_on_path_traced_back_from_node(
        collision_edge.detector_node, [&path_edges_1](const SearchGraphEdge& edge) {
            path_edges_1.push_back(edge);
        });

    std::vector<SearchGraphEdge> path_edges_2 = {collision_edge};
    auto other_node = collision_edge.detector_node->neighbors[collision_edge.neighbor_index];
    if (other_node)
        iter_edges_on_path_traced_back_from_node(other_node, [&path_edges_2](const SearchGraphEdge& edge) {
            path_edges_2.push_back(edge);
        });

    auto last_node_on_2 = path_edges_2.back().detector_node->neighbors[path_edges_2.back().neighbor_index];
    if (last_node_on_2 == loc_from_ptr) {
        // Reverse path 2
        reverse_path_and_handle_edges(path_edges_2, handle_edge);
        // Handle path 1 in the original order
        for (auto& x : path_edges_1)
            handle_edge(x);
    } else {
        // Reverse path 1
        reverse_path_and_handle_edges(path_edges_1, handle_edge);
        // Handle path 2 in the original order
        for (auto& x : path_edges_2)
            handle_edge(x);
    }

    reset();
}

}  // namespace pm

#endif  // PYMATCHING2_SEARCH_FLOODER_H
