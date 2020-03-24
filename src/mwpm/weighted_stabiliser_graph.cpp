#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include "weighted_stabiliser_graph.h"
#include <memory>
#include <stdexcept>


WeightedStabiliserGraph::WeightedStabiliserGraph(int num_stabilisers){
    wgraph_t sgraph = wgraph_t(num_stabilisers);
    this->stabiliser_graph = sgraph;
}

WeightedStabiliserGraph::WeightedStabiliserGraph(
            const py::array_t<int>& indices, 
            const py::array_t<double>& weights
){
    auto x = indices.unchecked<1>();
    assert((x.shape(0) % 2) == 0);
    int smax = 0;
    for (py::ssize_t i=0; i<(x.shape(0)); i++){
        if (x[i] > smax){
            smax = x[i];
        }
    }
    int num_stabilisers = smax+1;

    wgraph_t sgraph = wgraph_t(num_stabilisers);
    this->stabiliser_graph = sgraph;

    auto w = weights.unchecked<1>();
    
    assert(w.shape(0) == x.shape(0)/2);
    for (py::ssize_t i=0; i<x.shape(0)/2; i++){
        AddEdge(x[2*i], x[2*i+1], (int) i, w[i]);
    }
    ComputeAllPairsShortestPaths();
}

void WeightedStabiliserGraph::AddEdge(
    int node1, 
    int node2,
    int qubit_id,
    double weight){
        if (qubit_id < -1){
            throw std::runtime_error("Qubit ids must be non-negative, or -1 if the edge is not a qubit.");
        }
        boost::add_edge(
            boost::vertex(node1, stabiliser_graph), 
            boost::vertex(node2, stabiliser_graph), 
            {qubit_id, weight}, 
            stabiliser_graph);
}

void WeightedStabiliserGraph::ComputeAllPairsShortestPaths(){
    int n = boost::num_vertices(stabiliser_graph);
    all_distances.clear();
    all_predecessors.clear();
    for (int i=0; i<n; i++){
        std::vector<double> distances(n);
        std::vector<vertex_descriptor> p(n);
        vertex_descriptor from = boost::vertex(i, stabiliser_graph);
        boost::dijkstra_shortest_paths(stabiliser_graph, from,
            boost::weight_map(boost::get(&WeightedEdgeData::weight, stabiliser_graph))
            .distance_map(boost::make_iterator_property_map(distances.begin(),
                            boost::get(boost::vertex_index, stabiliser_graph)))
            .predecessor_map(&p[0]));
        all_distances.push_back(distances);
        all_predecessors.push_back(p);
    }
}

double WeightedStabiliserGraph::Distance(int node1, int node2) const {
    vertex_descriptor n2 = boost::vertex(node2, stabiliser_graph);
    return all_distances[node1][n2];
}

std::vector<int> WeightedStabiliserGraph::ShortestPath(int node1, int node2) const {
    std::vector<vertex_descriptor> parent = all_predecessors[node2];
    auto index = boost::get(boost::vertex_index, stabiliser_graph);
    int c = boost::vertex(node1, stabiliser_graph);
    std::vector<int> path;
    path.push_back(index[c]);
    while (parent[c]!=c){
        c = parent[c];
        path.push_back(index[c]);
    }
    return path;
}

int WeightedStabiliserGraph::GetNumQubits() const {
    auto qid = boost::get(&WeightedEdgeData::qubit_id, stabiliser_graph);
    int num_edges = boost::num_edges(stabiliser_graph);
    int num_qubits = 0;
    auto es = boost::edges(stabiliser_graph);
    for (auto eit = es.first; eit != es.second; ++eit) {
        int qubit = qid[*eit];
        if (qubit >=0){
            num_qubits++;
        } else if (qubit != -1){
            throw std::runtime_error("Qubit ids must be non-negative, or -1 if the edge is not a qubit.");
        }
    }
    return num_qubits;
}

int WeightedStabiliserGraph::GetNumStabilisers() const {
    return boost::num_vertices(stabiliser_graph);
};

int WeightedStabiliserGraph::QubitID(int node1, int node2) const {
    auto e = boost::edge(node1, node2, stabiliser_graph);
    if (!e.second){
        std::runtime_error("Graph does not contain edge (" 
                        + std::to_string((int)node1)
                        + std::to_string((int)node2) + ").");
    }
    return stabiliser_graph[e.first].qubit_id;
}