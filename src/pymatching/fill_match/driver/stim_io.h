#ifndef PYMATCHING2_STIM_IO_H
#define PYMATCHING2_STIM_IO_H

#include "pymatching/fill_match/flooder/graph.h"
#include "stim.h"

namespace pm {
pm::MatchingGraph detector_error_model_to_matching_graph(
    const stim::DetectorErrorModel &detector_error_model, pm::weight_int num_distinct_weights);

struct Neighbor {
    std::vector<Neighbor> *node;
    double probability;
    std::vector<size_t> observables;
};

class ProbabilityGraph {
   public:
    std::vector<std::vector<Neighbor>> nodes;
    size_t num_nodes;
    size_t num_observables;

    explicit ProbabilityGraph(size_t num_nodes, size_t num_observables)
        : num_nodes(num_nodes), num_observables(num_observables) {
        nodes.resize(num_nodes);
    };

    void add_or_merge_edge(size_t u, size_t v, double probability, const std::vector<size_t>& observables);

    void add_or_merge_boundary_edge(size_t u, double probability, const std::vector<size_t>& observables);

    void handle_dem_instruction(double p, const std::vector<size_t> &detectors, std::vector<size_t>& observables);

    pm::MatchingGraph to_matching_graph(pm::weight_int num_distinct_weights);

    double min_nonzero_probability();
};

}  // namespace pm

#endif  // PYMATCHING2_STIM_IO_H
