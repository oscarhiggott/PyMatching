#include "stim_io.h"
#include "graph.h"

#include <utility>
#include <algorithm>


void pm::ProbabilityGraph::add_or_merge_edge(size_t u, size_t v, double probability, pm::obs_int observables) {
    size_t larger_node = std::max(u, v);
    if (larger_node + 1 > nodes.size()) {
        throw std::invalid_argument("Node " + std::to_string(larger_node) + " exceeds number of nodes "
                                                                            "in graph (" + std::to_string(num_nodes) + ")");
    }
    auto it = std::find_if(nodes[u].begin(), nodes[u].end(),
                           [this, v](Neighbor& neighbor){
        return neighbor.node == &nodes[v];
    });
    if (it == nodes[u].end()){
        nodes[u].push_back({&nodes[v], probability, observables});
        nodes[v].push_back({&nodes[u], probability, observables});
    } else {
        double new_p = (*it).probability*(1-probability) + probability*(1-(*it).probability);
        (*it).probability = new_p;
        auto it2 = std::find_if((*it).node->begin(), (*it).node->end(),
                                [this, u](Neighbor& neighbor){
            return neighbor.node == &nodes[u];
        });
        if (it2 != (*it).node->end())
            (*it2).probability = new_p;
    }
}

void pm::ProbabilityGraph::add_or_merge_boundary_edge(size_t u, double probability, pm::obs_int observables) {
    if (u > nodes.size() - 1) {
        throw std::invalid_argument("Node " + std::to_string(u) + " exceeds number of nodes "
                                                                  "in graph (" + std::to_string(num_nodes) + ")");
    }
    auto it = std::find_if(nodes[u].begin(), nodes[u].end(),
                           [this](Neighbor& neighbor){
                               return neighbor.node == nullptr;
                           });
    if (it == nodes[u].end()){
        nodes[u].push_back({nullptr, probability, observables});
    } else {
        (*it).probability = (*it).probability*(1-probability) + probability*(1-(*it).probability);
    }
}

void
pm::ProbabilityGraph::handle_dem_instruction(double p, const std::vector<size_t> &detectors, pm::obs_int obs_mask) {
    if (detectors.size() == 2){
        add_or_merge_edge(detectors[0], detectors[1], p, obs_mask);
    } else if (detectors.size() == 1) {
        add_or_merge_boundary_edge(detectors[0], p, obs_mask);
    }
}

double pm::ProbabilityGraph::min_nonzero_probability() {
    double min_prob = 1.0;
    for (auto& node : nodes) {
        for (auto& neighbor : node) {
            if (neighbor.probability > 0 && neighbor.probability < min_prob)
                min_prob = neighbor.probability;
        }
    }
    return min_prob;
}

pm::MatchingGraph pm::ProbabilityGraph::to_matching_graph(pm::weight_int num_buckets) {
    double min_prob = min_nonzero_probability();
    double max_weight = std::log((1-min_prob)/min_prob);
    pm::MatchingGraph matching_graph(nodes.size());
    for (auto& node : nodes) {
        for (auto& neighbor : node) {
            auto i = &node - &nodes[0];
            double normed_weight = std::log((1-neighbor.probability)/neighbor.probability) / max_weight;
            pm::weight_int w = num_buckets * normed_weight;
            if (!neighbor.node){
                matching_graph.add_boundary_edge(i, w, neighbor.obs_mask);
            } else {
                auto j = neighbor.node - &nodes[0];
                matching_graph.add_edge(i, j, w, neighbor.obs_mask);
            }
        }
    }
    return matching_graph;
}


pm::MatchingGraph pm::detector_error_model_to_matching_graph(stim::DetectorErrorModel& detector_error_model,
                                                             pm::weight_int num_buckets) {
    pm::ProbabilityGraph probability_graph(detector_error_model.count_detectors());
    detector_error_model.iter_flatten_error_instructions(
            [&probability_graph](const stim::DemInstruction& instruction){
                std::vector<size_t> dets;
                pm::obs_int obs_mask = 0;
                double p = instruction.arg_data[0];
                for (auto& target : instruction.target_data) {
                    if (target.is_relative_detector_id()){
                        dets.push_back(target.val());
                    } else if (target.is_observable_id()) {
                        obs_mask ^= 1 << target.val();
                    } else if (target.is_separator()) {
                        probability_graph.handle_dem_instruction(p, dets, obs_mask);
                        obs_mask = 0;
                        dets.clear();
                    }
                }
                probability_graph.handle_dem_instruction(p, dets, obs_mask);
            }
            );

    return probability_graph.to_matching_graph(num_buckets);
}


