#pragma once

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <memory>
#include <set>
#include <cstdint>
#include "stabiliser_graph.h"

namespace py = pybind11;

struct WeightedEdgeData {
    std::set<int> qubit_ids;
    double weight;
    double error_probability;
    bool has_error_probability;
};

typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
    boost::no_property, WeightedEdgeData > wgraph_t;
typedef boost::graph_traits < wgraph_t >::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits < wgraph_t >::edge_descriptor edge_descriptor;


class WeightedStabiliserGraph : public IStabiliserGraph{
    public:
        wgraph_t stabiliser_graph;
        WeightedStabiliserGraph(
            int num_stabilisers,
            std::vector<int>& boundary
            );
        WeightedStabiliserGraph(
            const py::array_t<int>& indices, 
            const py::array_t<double>& weights,
            std::vector<int>& boundary
        );
        WeightedStabiliserGraph(
            const py::array_t<int>& indices, 
            const py::array_t<double>& weights,
            const py::array_t<double>& error_probabilies,
            std::vector<int>& boundary
        );
        void AddEdge(int node1, int node2, std::set<int> qubit_ids,
                     double weight, double error_probability,
                     bool has_error_probability);
        virtual void ComputeAllPairsShortestPaths();
        virtual double Distance(int node1, int node2);
        virtual std::vector<int> ShortestPath(int node1, int node2);
        virtual std::set<int> QubitIDs(int node1, int node2) const;
        virtual int GetNumQubits() const;
        virtual int GetNumStabilisers() const;
        std::pair<py::array_t<std::uint8_t>,py::array_t<std::uint8_t>> AddNoise() const;
        std::vector<std::vector<double>> all_distances;
        std::vector<std::vector<vertex_descriptor>> all_predecessors;
        bool all_edges_have_error_probabilities;
        virtual std::vector<int> GetBoundary() const;
        virtual void SetBoundary(std::vector<int>& boundary);
        virtual bool HasComputedAllPairsShortestPaths() const;
        void ResetDijkstraNeighbours();
        std::vector<std::pair<int, double>> GetNearestNeighbours(
            int source, int num_neighbours, std::vector<int>& defect_id);
        std::vector<int> GetPath(
            int source, int target);
        std::vector<int> boundary;
        std::vector<double> _distances;
        std::vector<int> _predecessors;
};