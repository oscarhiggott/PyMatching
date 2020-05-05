#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include "weighted_stabiliser_graph.h"
#include <memory>
#include <set>
#include <utility>
#include <stdexcept>
#include "rand_gen.h"


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
    double weight,
    double error_probability){
        if (qubit_id < -1){
            throw std::runtime_error("Qubit ids must be non-negative, or -1 if the edge is not a qubit.");
        }
        boost::add_edge(
            boost::vertex(node1, stabiliser_graph), 
            boost::vertex(node2, stabiliser_graph), 
            {qubit_id, weight, error_probability}, 
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
    int maxid = 0;
    std::set<int> qubits;
    auto es = boost::edges(stabiliser_graph);
    for (auto eit = es.first; eit != es.second; ++eit) {
        int qubit = qid[*eit];
        if (qubit >= maxid){
            maxid = qubit;
        }
        if (qubit >=0){
            qubits.insert(qubit);
        } else if (qubit != -1){
            throw std::runtime_error("Qubit ids must be non-negative, or -1 if the edge is not a qubit.");
        }
    }
    int num_qubits = qubits.size();
    if (maxid + 1 != num_qubits){
        throw std::runtime_error("Qubit ids must be numbered 0...(N-1).");
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

std::pair<py::array_t<int>,py::array_t<int>> WeightedStabiliserGraph::AddNoise() const {
    auto syndrome = new std::vector<int>(GetNumStabilisers(), 0);
    auto error = new std::vector<int>(GetNumQubits(), 0);
    double p;
    int qid;
    vertex_descriptor s, t;
    bool to_flip;
    auto es = boost::edges(stabiliser_graph);
    for (auto eit = es.first; eit != es.second; ++eit) {
        p = stabiliser_graph[*eit].error_probability;
        if ((p >= 0) && (rand_float(0.0, 1.0) < p)){
            s = boost::source(*eit, stabiliser_graph);
            t = boost::target(*eit, stabiliser_graph);
            (*syndrome)[s] = ((*syndrome)[s] + 1) % 2;
            (*syndrome)[t] = ((*syndrome)[t] + 1) % 2;
            qid = stabiliser_graph[*eit].qubit_id;
            if (qid >= 0){
                (*error)[qid] = ((*error)[qid] + 1) % 2;
            }
        }
    }
    auto capsule = py::capsule(syndrome, [](void *syndrome) { delete reinterpret_cast<std::vector<int>*>(syndrome); });
    py::array_t<int> syndrome_arr = py::array_t<int>(syndrome->size(), syndrome->data(), capsule);
    auto err_capsule = py::capsule(error, [](void *error) { delete reinterpret_cast<std::vector<int>*>(error); });
    py::array_t<int> error_arr = py::array_t<int>(error->size(), error->data(), err_capsule);
    return {error_arr, syndrome_arr};
}