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
            const py::array_t<double>& weights, 
            int num_stabilisers, 
            int num_qubits
){
    wgraph_t sgraph = wgraph_t(num_stabilisers);
    this->stabiliser_graph = sgraph;
    auto x = indices.unchecked<1>();
    auto w = weights.unchecked<1>();
    assert((x.shape(0) % 2) ==0 );
    assert(x.shape(0)/2 == num_qubits);
    assert(w.shape(0) == num_qubits);
    for (py::ssize_t i=0; i<num_qubits; i++){
        AddEdge(x[2*i], x[2*i+1], (int) i, w[i]);
    }
    ComputeAllPairsShortestPaths();
}

void WeightedStabiliserGraph::AddEdge(
    int node1, 
    int node2,
    int qubit_id,
    double weight){
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
            boost::weight_map(get(&WeightedEdgeData::weight, stabiliser_graph))
            .distance_map(boost::make_iterator_property_map(distances.begin(),
                            boost::get(boost::vertex_index, stabiliser_graph)))
            .predecessor_map(&p[0]));
        all_distances.push_back(distances);
        all_predecessors.push_back(p);
    }
}

int WeightedStabiliserGraph::Distance(int node1, int node2) const {
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
    return boost::num_edges(stabiliser_graph);
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