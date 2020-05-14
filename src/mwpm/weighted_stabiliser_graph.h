#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <memory>
#include "stabiliser_graph.h"

namespace py = pybind11;

struct WeightedEdgeData {
    int qubit_id;
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
            int boundary = -1
            );
        WeightedStabiliserGraph(
            const py::array_t<int>& indices, 
            const py::array_t<double>& weights,
            int boundary = -1
        );
        WeightedStabiliserGraph(
            const py::array_t<int>& indices, 
            const py::array_t<double>& weights,
            const py::array_t<double>& error_probabilies,
            int boundary = -1
        );
        void AddEdge(int node1, int node2, int qubit_id,
                     double weight, double error_probability,
                     bool has_error_probability);
        void ComputeAllPairsShortestPaths();
        virtual double Distance(int node1, int node2) const;
        virtual std::vector<int> ShortestPath(int node1, int node2) const;
        virtual int QubitID(int node1, int node2) const;
        virtual int GetNumQubits() const;
        virtual int GetNumStabilisers() const;
        std::pair<py::array_t<int>,py::array_t<int>> AddNoise() const;
        std::vector<std::vector<double>> all_distances;
        std::vector<std::vector<vertex_descriptor>> all_predecessors;
        bool all_edges_have_error_probabilities;
        virtual int GetBoundaryVertex() const;
        virtual void SetBoundaryVertex(int boundary);
        int boundary;
};