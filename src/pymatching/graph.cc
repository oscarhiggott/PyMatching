#include "pymatching/graph.h"


#include "pymatching/graph_fill_region.h"
#include "pymatching/events.h"

namespace pm {

void MatchingGraph::add_edge(size_t u, size_t v, weight_int weight, obs_int observables) {
    size_t larger_node = std::max(u, v);
    if (larger_node + 1 > nodes.size()) {
        throw std::invalid_argument(
            "Node " + std::to_string(larger_node) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }

    // Allow parallel edges?
    nodes[u].neighbors.push_back(&(nodes[v]));
    nodes[u].neighbor_weights.push_back(weight);
    nodes[u].neighbor_observables.push_back(observables);

    nodes[v].neighbors.push_back(&(nodes[u]));
    nodes[v].neighbor_weights.push_back(weight);
    nodes[v].neighbor_observables.push_back(observables);
}

void MatchingGraph::add_boundary_edge(size_t u, weight_int weight, obs_int observables) {
    if (u >= nodes.size()) {
        throw std::invalid_argument(
            "Node " + std::to_string(u) +
            " exceeds number of nodes "
            "in graph (" +
            std::to_string(num_nodes) + ")");
    }
    auto &n = nodes[u];
    if (!n.neighbors.empty() && n.neighbors[0] == nullptr) {
        throw std::invalid_argument("Max one boundary edge.");
    }
    n.neighbors.insert(n.neighbors.begin(), 1, nullptr);
    n.neighbor_weights.insert(n.neighbor_weights.begin(), 1, weight);
    n.neighbor_observables.insert(n.neighbor_observables.begin(), 1, observables);
}

MatchingGraph::MatchingGraph(size_t num_nodes) : num_nodes(num_nodes) {
    nodes.resize(num_nodes);
}

MatchingGraph::MatchingGraph(MatchingGraph &&graph) noexcept
    : nodes(std::move(graph.nodes)), num_nodes(graph.num_nodes) {
}

MatchingGraph::MatchingGraph() : num_nodes(0) {
}

Varying32 DetectorNode::total_radius() const {
    if (!reached_from_source)
        return pm::Varying32(0);
    auto curr_region = reached_from_source->region_that_arrived;
    decltype(curr_region->radius.data) tot_rad = 0;
    while (curr_region->blossom_parent) {
        tot_rad += curr_region->radius.y_intercept();
        curr_region = curr_region->blossom_parent;
    }
    return curr_region->radius + tot_rad;
}

Varying32 DetectorNode::local_radius() const {
    if (!reached_from_source)
        return Varying32(0);
    return total_radius() - distance_from_source;
}

bool DetectorNode::has_same_owner_as(const DetectorNode &other) const {
    if ((region_that_arrived == nullptr) != (other.region_that_arrived == nullptr))
        return false;
    if (region_that_arrived == other.region_that_arrived)
        return true;
    return top_region() == other.top_region();
}

GraphFillRegion *DetectorNode::top_region() const {
    if (!region_that_arrived)
        return nullptr;
    return region_that_arrived->top_region();
}

void DetectorNode::reset() {
    observables_crossed_from_source = 0;
    reached_from_source = nullptr;
    distance_from_source = 0;
    region_that_arrived = nullptr;
    node_event_tracker.clear();
}

size_t DetectorNode::index_of_neighbor(DetectorNode *target) const {
    for (size_t k = 0; k < neighbors.size(); k++) {
        if (neighbors[k] == target) {
            return k;
        }
    }
    throw std::invalid_argument("Failed to find neighbor.");
}

}  // namespace pm
