#ifndef PYMATCHING2_GRAPH_H
#define PYMATCHING2_GRAPH_H


#include <vector>
#include "fixed_length_vector.h"
#include "varying.h"


namespace pm {

    typedef uint64_t obs_int;
    typedef uint8_t weight_int;

    class GraphFillRegion;

    struct TentativeEvent;

//    template<template<class> class neighbour_list>
    class DetectorNode{
    public:
        DetectorNode() :
            observables_crossed_from_source(0),
            reached_from_source(nullptr),
            distance_from_source(0),
            region_that_arrived(nullptr) {}


        obs_int observables_crossed_from_source;
        DetectorNode* reached_from_source;
        time_int distance_from_source;
        GraphFillRegion* region_that_arrived;

        std::vector<DetectorNode*> neighbors;
        std::vector<weight_int> neighbor_weights;
        std::vector<obs_int> neighbor_observables;
        std::vector<TentativeEvent*> neighbor_schedules;
    };


    class Graph {
    public:
        size_t num_nodes;
        Graph() = delete;
        explicit Graph(size_t num_nodes);
        std::vector<DetectorNode> nodes;
        void add_edge(size_t u, size_t v, weight_int weight, obs_int observables);
        void add_boundary_edge(size_t u, weight_int weight, obs_int observables);
    };

}


#endif //PYMATCHING2_GRAPH_H
