// Copyright 2020 Oscar Higgott

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//      http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include <boost/graph/connected_components.hpp>
#include "weighted_stabiliser_graph.h"
#include <memory>
#include <set>
#include <utility>
#include <stdexcept>
#include <cstdint>
#include <iostream>
#include <limits>
#include "rand_gen.h"


WeightedMatchingGraph::WeightedMatchingGraph()
    : all_edges_have_error_probabilities(true),
     connected_components_need_updating(true) {
    wgraph_t sgraph = wgraph_t();
    this->stabiliser_graph = sgraph;
}


WeightedMatchingGraph::WeightedMatchingGraph(
    int num_stabilisers,
    std::set<int>& boundary)
    : all_edges_have_error_probabilities(true),
    boundary(boundary),
    connected_components_need_updating(true) {
    wgraph_t sgraph = wgraph_t(num_stabilisers+boundary.size());
    this->stabiliser_graph = sgraph;
}

void WeightedMatchingGraph::AddEdge(
    int node1, 
    int node2, 
    std::set<int> qubit_ids, 
    double weight, 
    double error_probability, 
    bool has_error_probability){
        if (!has_error_probability){
            all_edges_have_error_probabilities = false;
        }
        connected_components_need_updating = true;
        WeightedEdgeData data;
        data.qubit_ids = qubit_ids;
        data.weight = weight;
        data.error_probability = error_probability;
        data.has_error_probability = has_error_probability;
        boost::add_edge(
            boost::vertex(node1, stabiliser_graph), 
            boost::vertex(node2, stabiliser_graph), 
            data, 
            stabiliser_graph);
}

void WeightedMatchingGraph::ComputeAllPairsShortestPaths(){
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

class exit_search{};


class DijkstraNeighbourVisitor : public boost::default_dijkstra_visitor
{
    public:
        DijkstraNeighbourVisitor(std::vector<int>& defect_id, 
            int num_defects, std::vector<int>& examined_defects,
            std::vector<int>& discovered_nodes) : 
        defect_id(defect_id), num_defects(num_defects), 
        examined_defects(examined_defects), num_found(0),
        discovered_nodes(discovered_nodes) {}

        void examine_vertex(wgraph_t::vertex_descriptor v, const wgraph_t &g)
        {
            if (defect_id[v] > -1){
                num_found++;
                examined_defects.push_back(v);
                if (num_found >= num_defects) {
                    throw exit_search();
                }
            }   
        }
        void discover_vertex(wgraph_t::vertex_descriptor v, const wgraph_t &g){
            discovered_nodes.push_back(v);
        }
        std::vector<int>& defect_id;
        int num_defects;
        std::vector<int>& examined_defects;
        std::vector<int>& discovered_nodes;
        int num_found;
};


void WeightedMatchingGraph::ResetDijkstraNeighbours(){
    int n = boost::num_vertices(stabiliser_graph);
    double inf = std::numeric_limits<double>::max();
    if (_distances.size() < n){
        _distances.resize(n, inf);
    }
    if (_predecessors.size() < n){
        _predecessors.resize(n);
        for (int i=0; i<_predecessors.size(); i++){
            _predecessors[i] = i;
        }
    }
}


std::vector<std::pair<int, double>> WeightedMatchingGraph::GetNearestNeighbours(
    int source, int num_neighbours, std::vector<int>& defect_id){
    int n = boost::num_vertices(stabiliser_graph);
    assert(source < n);
    assert(defect_id.size() == n);
    double inf = std::numeric_limits<double>::max();
    ResetDijkstraNeighbours();
    _distances[source] = 0;
    std::vector<int> examined_defects;
    std::vector<int> discovered_nodes;
    DijkstraNeighbourVisitor vis = DijkstraNeighbourVisitor(
        defect_id, num_neighbours, examined_defects, discovered_nodes);
    vertex_descriptor from = boost::vertex(source, stabiliser_graph);
    try {
        boost::dijkstra_shortest_paths_no_color_map_no_init(stabiliser_graph, from, 
                &_predecessors[0], boost::make_iterator_property_map(_distances.begin(),
                boost::get(boost::vertex_index, stabiliser_graph)), 
                boost::get(&WeightedEdgeData::weight, stabiliser_graph),
                boost::get(boost::vertex_index, stabiliser_graph), 
                std::less<double>(),
                boost::closed_plus<double>(),
                inf,
                0,
                vis);
    } catch (exit_search e) {}
    
    std::vector<std::pair<int, double>> neighbours;
    for (auto d : examined_defects){
        neighbours.push_back({d, _distances[d]});
    }

    for (auto n : discovered_nodes){
        _distances[n] = inf;
        _predecessors[n] = n;
    }
    return neighbours;    
}


class DijkstraPathVisitor : public boost::default_dijkstra_visitor
{
    public:
        DijkstraPathVisitor(int target, std::vector<int>& discovered_nodes) : 
        target(target),
        discovered_nodes(discovered_nodes) {}

        void examine_vertex(wgraph_t::vertex_descriptor v, const wgraph_t &g)
        {
            if (v == target){
                throw exit_search();
            }   
        }
        void discover_vertex(wgraph_t::vertex_descriptor v, const wgraph_t &g){
            discovered_nodes.push_back(v);
        }
        int target;
        std::vector<int>& discovered_nodes;
};


std::vector<int> WeightedMatchingGraph::GetPath(
    int source, int target){
    int n = boost::num_vertices(stabiliser_graph);
    assert(source < n);
    assert(target < n);
    double inf = std::numeric_limits<double>::max();
    ResetDijkstraNeighbours();
    _distances[source] = 0;
    std::vector<int> discovered_nodes;
    DijkstraPathVisitor vis = DijkstraPathVisitor(
        target, discovered_nodes);
    vertex_descriptor from = boost::vertex(source, stabiliser_graph);
    try {
        boost::dijkstra_shortest_paths_no_color_map_no_init(stabiliser_graph, from, 
                &_predecessors[0], boost::make_iterator_property_map(_distances.begin(),
                boost::get(boost::vertex_index, stabiliser_graph)), 
                boost::get(&WeightedEdgeData::weight, stabiliser_graph),
                boost::get(boost::vertex_index, stabiliser_graph), 
                std::less<double>(),
                boost::closed_plus<double>(),
                inf,
                0,
                vis);
    } catch (exit_search e) {}

    std::vector<int> path;
    path.push_back(target);
    while (_predecessors[target] != target){
        target = _predecessors[target];
        path.push_back(target);
    }
    std::reverse(path.begin(),path.end());

    for (auto n : discovered_nodes){
        _distances[n] = inf;
        _predecessors[n] = n;
    }
    return path;    
}


double WeightedMatchingGraph::Distance(int node1, int node2) {
    if (!HasComputedAllPairsShortestPaths()){
        ComputeAllPairsShortestPaths();
    }
    vertex_descriptor n2 = boost::vertex(node2, stabiliser_graph);
    return all_distances[node1][n2];
}

std::vector<int> WeightedMatchingGraph::ShortestPath(int node1, int node2) {
    if (!HasComputedAllPairsShortestPaths()){
        ComputeAllPairsShortestPaths();
    }
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


int WeightedMatchingGraph::GetNumEdges() const {
    return boost::num_edges(stabiliser_graph);
}


int WeightedMatchingGraph::GetNumQubits() const {
    auto qid = boost::get(&WeightedEdgeData::qubit_ids, stabiliser_graph);
    int num_edges = boost::num_edges(stabiliser_graph);
    int maxid = -1;
    std::set<int> qubits;
    std::set<int> edge_qubits;
    auto es = boost::edges(stabiliser_graph);
    for (auto eit = es.first; eit != es.second; ++eit) {
        edge_qubits = qid[*eit];
        for (auto qubit : edge_qubits){
            if (qubit >= maxid){
                maxid = qubit;
            }
            if (qubit >=0){
                qubits.insert(qubit);
            } else if (qubit != -1){
                throw std::runtime_error("Qubit ids must be non-negative, or -1 if the edge is not a qubit.");
            }
        }
    }
    int num_qubits = qubits.size();
    if (maxid + 1 != num_qubits){
        throw std::runtime_error("Qubit ids must be numbered 0...(N-1).");
    }
    return num_qubits;
}

int WeightedMatchingGraph::GetNumNodes() const {
    return boost::num_vertices(stabiliser_graph);
};

std::set<int> WeightedMatchingGraph::QubitIDs(int node1, int node2) const {
    auto e = boost::edge(node1, node2, stabiliser_graph);
    if (!e.second){
        std::runtime_error("Graph does not contain edge (" 
                        + std::to_string((int)node1)
                        + std::to_string((int)node2) + ").");
    }
    return stabiliser_graph[e.first].qubit_ids;
}

std::pair<py::array_t<std::uint8_t>,py::array_t<std::uint8_t>> WeightedMatchingGraph::AddNoise() const {
    auto syndrome = new std::vector<int>(GetNumNodes(), 0);
    auto error = new std::vector<int>(GetNumQubits(), 0);
    double p;
    std::set<int> qids;
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
            qids = stabiliser_graph[*eit].qubit_ids;
            for (auto qid : qids){
                if (qid >= 0){
                    (*error)[qid] = ((*error)[qid] + 1) % 2;
                }
            }
        }
    }
    for (auto b : boundary){
        (*syndrome)[b] = 0;
    }

    auto capsule = py::capsule(syndrome, [](void *syndrome) { delete reinterpret_cast<std::vector<int>*>(syndrome); });
    py::array_t<int> syndrome_arr = py::array_t<int>(syndrome->size(), syndrome->data(), capsule);
    auto err_capsule = py::capsule(error, [](void *error) { delete reinterpret_cast<std::vector<int>*>(error); });
    py::array_t<int> error_arr = py::array_t<int>(error->size(), error->data(), err_capsule);
    return {error_arr, syndrome_arr};
}

std::set<int> WeightedMatchingGraph::GetBoundary() const {
    return boundary;
}

void WeightedMatchingGraph::SetBoundary(std::set<int>& boundary) {
    this->boundary = boundary;
    connected_components_need_updating = true;
    return;
}

std::vector<std::tuple<int,int,WeightedEdgeData>> WeightedMatchingGraph::GetEdges() const {
    std::vector<std::tuple<int,int,WeightedEdgeData>> edges;
    auto es = boost::edges(stabiliser_graph);
    for (auto eit = es.first; eit != es.second; ++eit) {
        WeightedEdgeData edata = stabiliser_graph[*eit];
        int s = boost::source(*eit, stabiliser_graph);
        int t = boost::target(*eit, stabiliser_graph);
        std::tuple<int,int,WeightedEdgeData> edge = std::make_tuple(s, t, edata);
        edges.push_back(edge);
    }
    return edges;
}

bool WeightedMatchingGraph::HasComputedAllPairsShortestPaths() const {
    int n = boost::num_vertices(stabiliser_graph);
    bool has_distances = all_distances.size() == n;
    bool has_preds = all_predecessors.size() == n;
    return has_distances && has_preds;
}

int WeightedMatchingGraph::GetNumConnectedComponents() {
    if (connected_components_need_updating){
        component.resize(boost::num_vertices(stabiliser_graph));
        num_components = boost::connected_components(stabiliser_graph, &component[0]);
        component_boundary.resize(num_components, -1);
        for (auto b : boundary){
            int c = component[b];
            if (component_boundary[c] == -1){
                component_boundary[c] = b;
            }
        }
    }
    connected_components_need_updating = false;
    return num_components;
}

void WeightedMatchingGraph::FlipBoundaryNodesIfNeeded(std::set<int> &defects){
    int num_comps = GetNumConnectedComponents();
    if (num_comps == 1){
        if ((defects.size() % 2) == 0){
            return;
        } else if ((defects.size() % 2) == 1 && boundary.size() == 0){
            throw std::invalid_argument(
            "The syndrome has an odd number of defects, but no boundary nodes were provided"
            );
        }
    }
    std::vector<std::uint8_t> component_parities(num_comps, 0);
    for (auto df : defects){
        if (df >= component.size()){
            throw std::invalid_argument(
            "Defect id should not exceed the number of vertices in the graph"
            );
        }
        component_parities[component[df]] ^= 1;
    }

    for (int i=0; i<component_parities.size(); i++){
        if (component_parities[i] == 1){
            int b = component_boundary[i];
            if (b == -1){
                throw std::invalid_argument(
                    "The syndrome has an odd number of defects in a component of the matching graph "
                    "that does not have a boundary node"
                    );
            }
            bool is_in_defects = defects.find(b) != defects.end();
            if (is_in_defects){
                defects.erase(b);
            } else {
                defects.insert(b);
            }
        }
    }
}

bool WeightedMatchingGraph::AllEdgesHaveErrorProbabilities() const {
    return all_edges_have_error_probabilities;
}