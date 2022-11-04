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

#include "search_flooder.h"

#include <limits>

pm::SearchFlooder::SearchFlooder() : target_type(NO_TARGET) {
}

pm::SearchFlooder::SearchFlooder(pm::SearchGraph graph) : graph(std::move(graph)), target_type(NO_TARGET) {
}

std::pair<size_t, pm::cumulative_time_int> pm::SearchFlooder::find_next_event_at_node_returning_neighbor_index_and_time(
    const pm::SearchDetectorNode &detector_node) const {
    pm::cumulative_time_int best_time = std::numeric_limits<cumulative_time_int>::max();
    size_t best_neighbor = SIZE_MAX;

    size_t start = 0;
    if (!detector_node.neighbors.empty() && detector_node.neighbors[0] == nullptr) {
        if (target_type == BOUNDARY) {
            auto weight = detector_node.neighbor_weights[0];
            auto collision_time = queue.cur_time + weight - (queue.cur_time - detector_node.distance_from_source);
            if (collision_time < best_time) {
                best_time = collision_time;
                best_neighbor = 0;
            }
        }
        start++;
    }

    // Handle non-boundary neighbors.
    for (size_t i = start; i < detector_node.neighbors.size(); i++) {
        auto weight = detector_node.neighbor_weights[i];
        auto neighbor = detector_node.neighbors[i];

        pm::cumulative_time_int collision_time;

        if (neighbor->reached_from_source == detector_node.reached_from_source) {
            continue;
        } else if (!neighbor->reached_from_source) {
            collision_time = queue.cur_time + weight - (queue.cur_time - detector_node.distance_from_source);
        } else {
            // neighbor->reached_from_source != detector_node.reached_from_source
            auto covered_from_this_node = queue.cur_time - detector_node.distance_from_source;
            auto covered_from_neighbor = queue.cur_time - neighbor->distance_from_source;
            collision_time = queue.cur_time + (weight - covered_from_this_node - covered_from_neighbor) / 2;
        }

        if (collision_time < best_time) {
            best_time = collision_time;
            best_neighbor = i;
        }
    }

    return {best_neighbor, best_time};
}

void pm::SearchFlooder::reschedule_events_at_search_detector_node(pm::SearchDetectorNode &detector_node) {
    auto x = find_next_event_at_node_returning_neighbor_index_and_time(detector_node);
    if (x.first == SIZE_MAX) {
        detector_node.node_event_tracker.set_no_desired_event();
    } else {
        detector_node.node_event_tracker.set_desired_event({&detector_node, cyclic_time_int{x.second}}, queue);
    }
}

void pm::SearchFlooder::do_search_starting_at_empty_search_detector_node(pm::SearchDetectorNode *src) {
    src->reached_from_source = src;
    src->index_of_predecessor = SIZE_MAX;
    src->distance_from_source = 0;
    reached_nodes.push_back(src);
    reschedule_events_at_search_detector_node(*src);
}

void pm::SearchFlooder::do_search_exploring_empty_detector_node(
    pm::SearchDetectorNode &empty_node, size_t empty_to_from_index) {
    auto from_node = empty_node.neighbors[empty_to_from_index];
    empty_node.reached_from_source = from_node->reached_from_source;
    empty_node.index_of_predecessor = empty_to_from_index;
    empty_node.distance_from_source =
        empty_node.neighbor_weights[empty_to_from_index] + from_node->distance_from_source;
    reached_nodes.push_back(&empty_node);
    reschedule_events_at_search_detector_node(empty_node);
}

pm::SearchGraphEdge pm::SearchFlooder::do_look_at_node_event(pm::SearchDetectorNode &node) {
    auto next = find_next_event_at_node_returning_neighbor_index_and_time(node);
    if (next.second == queue.cur_time) {
        auto dst = node.neighbors[next.first];
        if (!dst) {
            // Found boundary (and must have target_type == BOUNDARY). Return collision edge for search to terminate.
            return {&node, next.first};
        } else if (node.reached_from_source && !dst->reached_from_source) {
            do_search_exploring_empty_detector_node(*dst, dst->index_of_neighbor(&node));
            // Need to revisit node immediately, since the search hasn't finished and other events may need to be
            // handled along other adjacent edges.
            node.node_event_tracker.set_desired_event(
                {
                    &node,
                    cyclic_time_int{queue.cur_time},
                },
                queue);
            return {nullptr, SIZE_MAX};
        } else {
            return {&node, next.first};
        }
    } else if (next.first != SIZE_MAX) {
        // Need to revisit this node, but at a later time
        node.node_event_tracker.set_desired_event(
            {
                &node,
                cyclic_time_int{next.second},
            },
            queue);
    }
    return {nullptr, SIZE_MAX};
}

pm::SearchGraphEdge pm::SearchFlooder::run_until_collision(pm::SearchDetectorNode *src, pm::SearchDetectorNode *dst) {
    do_search_starting_at_empty_search_detector_node(src);
    if (!dst) {
        target_type = BOUNDARY;
    } else {
        do_search_starting_at_empty_search_detector_node(dst);
        target_type = DETECTOR_NODE;
    }

    while (!queue.empty()) {
        FloodCheckEvent ev = queue.dequeue();
        if (ev.data_look_at_search_node->node_event_tracker.dequeue_decision(ev, queue)) {
            auto collision_edge = do_look_at_node_event(*ev.data_look_at_search_node);
            if (collision_edge.detector_node) {
                return collision_edge;
            }
        }
    }
    return {nullptr, SIZE_MAX};
}

void pm::SearchFlooder::reset_graph() {
    for (auto &detector_node : reached_nodes)
        detector_node->reset();
    reached_nodes.clear();
}

void pm::SearchFlooder::reset() {
    reset_graph();
    queue.reset();
}

pm::SearchFlooder::SearchFlooder(pm::SearchFlooder &&other) noexcept
    : graph(std::move(other.graph)),
      queue(std::move(other.queue)),
      reached_nodes(std::move(other.reached_nodes)),
      target_type(other.target_type) {
}
