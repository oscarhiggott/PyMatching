#ifndef PYMATCHING2_STIM_IO_H
#define PYMATCHING2_STIM_IO_H

#include "graph.h"
#include "stim.h"

namespace pm {
pm::MatchingGraph detector_error_model_to_matching_graph(
    const stim::DetectorErrorModel &detector_error_model, pm::weight_int num_buckets);

struct Neighbor {
    std::vector<Neighbor> *node;
    double probability;
    pm::obs_int obs_mask;
};

class ProbabilityGraph {
   public:
    std::vector<std::vector<Neighbor>> nodes;
    size_t num_nodes;

    explicit ProbabilityGraph(size_t num_nodes) : num_nodes(num_nodes) {
        nodes.resize(num_nodes);
    };

    void add_or_merge_edge(size_t u, size_t v, double probability, pm::obs_int observables);

    void add_or_merge_boundary_edge(size_t u, double probability, pm::obs_int observables);

    void handle_dem_instruction(double p, const std::vector<size_t> &detectors, pm::obs_int obs_mask);

    pm::MatchingGraph to_matching_graph(pm::weight_int num_buckets);

    double min_nonzero_probability();
};

}  // namespace pm

#endif  // PYMATCHING2_STIM_IO_H
