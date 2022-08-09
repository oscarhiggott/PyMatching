#ifndef PYMATCHING2_STIM_IO_H
#define PYMATCHING2_STIM_IO_H

#include "stim.h"
#include "graph.h"


//stim::DetectorErrorModel::from_file();
//auto reader = stim::MeasureRecordReader::make(nullptr, /*stim::format_name_to_enum_map*/ stim::SAMPLE_FORMAT_B8, 0 ,0, 0);
//stim::SparseShot sparse_shot;
//reader->start_and_read_entire_record(sparse_shot);
////        stim::find_argument()

namespace pm {
    pm::MatchingGraph detector_error_model_to_matching_graph(stim::DetectorErrorModel &detector_error_model,
                                                             pm::weight_int num_buckets = 1000);

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

        void handle_dem_instruction(double p, const std::vector<size_t>& detectors, pm::obs_int obs_mask);

        pm::MatchingGraph to_matching_graph(pm::weight_int num_buckets);

        double min_nonzero_probability();
    };


}



#endif //PYMATCHING2_STIM_IO_H
