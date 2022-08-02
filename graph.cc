#include<algorithm>
#include "graph.h"
#include "events.h"

namespace pm {

    void Graph::add_edge(size_t u, size_t v, weight_int weight, obs_int observables) {
        size_t larger_node = std::max(u, v);
        if (larger_node + 1 > nodes.size()) {
            throw std::invalid_argument("Node " + std::to_string(larger_node) + " exceeds number of nodes "
                                "in graph (" + std::to_string(num_nodes) + ")");
        }

        // Allow parallel edges?
        nodes[u].neighbors.push_back(&(nodes[v]));
        nodes[u].neighbor_weights.push_back(weight);
        nodes[u].neighbor_observables.push_back(observables);
        nodes[u].neighbor_schedules.push_back(nullptr);

        nodes[v].neighbors.push_back(&(nodes[u]));
        nodes[v].neighbor_weights.push_back(weight);
        nodes[v].neighbor_observables.push_back(observables);
        nodes[v].neighbor_schedules.push_back(nullptr);
    }

    void Graph::add_boundary_edge(size_t u, weight_int weight, obs_int observables) {
        if (u > nodes.size() - 1) {
            throw std::invalid_argument("Node " + std::to_string(u) + " exceeds number of nodes "
                          "in graph (" + std::to_string(num_nodes) + ")");
        }
        nodes[u].neighbors.push_back(nullptr);
        nodes[u].neighbor_weights.push_back(weight);
        nodes[u].neighbor_observables.push_back(observables);
        nodes[u].neighbor_schedules.push_back(nullptr);
    }

    Graph::Graph(size_t num_nodes) : num_nodes(num_nodes) {
        nodes.resize(num_nodes);
    }

    Graph::Graph(Graph &&graph)  noexcept : nodes(std::move(graph.nodes)), num_nodes(graph.num_nodes) {}

    void DetectorNode::invalidate_involved_schedule_items() {
        for (auto & neighbor_schedule : neighbor_schedules) {
            if (neighbor_schedule)
                neighbor_schedule->invalidate();
        }
    }

    Varying32 DetectorNode::total_radius() const {
        if (!reached_from_source)
            return pm::Varying32(0);
        auto curr_region = reached_from_source->region_that_arrived;
        typeof(curr_region->radius.data) tot_rad = 0;
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
        if ((!region_that_arrived) != (!other.region_that_arrived))
            return false;
        if (region_that_arrived == other.region_that_arrived)
            return true;
        return top_region() == other.top_region();
    }

    GraphFillRegion * DetectorNode::top_region() const {
        if (!region_that_arrived)
            return nullptr;
        return region_that_arrived->top_region();
    }

}
