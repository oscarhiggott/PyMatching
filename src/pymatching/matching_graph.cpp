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
#include "matching_graph.h"
#include <memory>
#include <set>
#include <utility>
#include <stdexcept>
#include <cstdint>
#include <iostream>
#include <limits>
#include <cmath>
#include "rand_gen.h"


WeightedEdgeData::WeightedEdgeData() {}

WeightedEdgeData::WeightedEdgeData(
    std::set<int> fault_ids,
    double weight,
    double error_probability,
    bool has_error_probability,
    bool weight_is_negative
): fault_ids(fault_ids), weight(weight),
error_probability(error_probability), has_error_probability(has_error_probability),
weight_is_negative(weight_is_negative) {}


std::string set_repr(std::set<int> x) {
    std::stringstream ss;
    ss << "{";
    bool first = true;
    for (auto i : x){
        if (first){
            first = false;
        } else {
            ss << ", ";
        }
        ss << i;
    }
    ss << "}";
    return ss.str();
}


std::string WeightedEdgeData::repr() const {
    std::stringstream ss;
    ss << "pymatching._cpp_mwpm.WeightedEdgeData(";
    ss << set_repr(fault_ids) << ", " << weight << ", " << error_probability << ", "
    << has_error_probability << ")";
    return ss.str();
}


MatchingGraph::MatchingGraph()
    : all_edges_have_error_probabilities(true),
     connected_components_need_updating(true),
     negative_weight_sum(0.0)  {
    wgraph_t sgraph = wgraph_t();
    this->matching_graph = sgraph;
}


MatchingGraph::MatchingGraph(
    int num_detectors,
    std::set<int>& boundary)
    : all_edges_have_error_probabilities(true),
    boundary(boundary),
    connected_components_need_updating(true),
     negative_weight_sum(0.0) {
    wgraph_t sgraph = wgraph_t(num_detectors+boundary.size());
    this->matching_graph = sgraph;
}

void MatchingGraph::AddEdge(
    int node1, 
    int node2, 
    std::set<int> fault_ids,
    double weight, 
    double error_probability, 
    bool has_error_probability){
        if (has_error_probability && (error_probability > 1 || error_probability < 0)){
            throw std::invalid_argument("error_probability must be between 0 and 1");
        }
        if (node1 < 0 || node2 < 0){
            throw std::invalid_argument("Node IDs must be non-negative");
        }
        auto n1 = boost::vertex(node1, matching_graph);
        auto n2 = boost::vertex(node2, matching_graph);
        int num_nodes = GetNumNodes();
        bool nodes_in_graph = (n1 < num_nodes) && (n2 < num_nodes);
        if (nodes_in_graph && boost::edge(n1, n2, matching_graph).second){
            throw std::invalid_argument("This edge already exists in the graph. "
                                        "Parallel edges are not supported.");
        }
        if (std::signbit(weight)){
            HandleNewNegativeWeightEdge(node1, node2, weight, fault_ids);
        }
        if (!has_error_probability){
            all_edges_have_error_probabilities = false;
            error_probability = -1;
        }
        connected_components_need_updating = true;
        WeightedEdgeData data;
        data.fault_ids = fault_ids;
        data.weight = std::abs(weight);
        data.error_probability = error_probability;
        data.has_error_probability = has_error_probability;
        data.weight_is_negative = std::signbit(weight);
        boost::add_edge(
            n1,
            n2,
            data, 
            matching_graph);
}


void MatchingGraph::HandleNewNegativeWeightEdge(int u, int v, double weight, std::set<int> &fault_ids){
    assert(std::signbit(weight));
    negative_weight_sum += weight;

    for (auto fid : fault_ids){
        if (negative_edge_fault_ids.find(fid) != negative_edge_fault_ids.end()){
            negative_edge_fault_ids.erase(fid);
        } else {
            negative_edge_fault_ids.insert(fid);
        }
    }

    for (auto node : {u, v}){
        if (negative_edge_syndrome.find(node) != negative_edge_syndrome.end()){
            negative_edge_syndrome.erase(node);
        } else {
            negative_edge_syndrome.insert(node);
        }
    }

}


void MatchingGraph::ComputeAllPairsShortestPaths(){
    int n = boost::num_vertices(matching_graph);
    all_distances.clear();
    all_predecessors.clear();
    for (int i=0; i<n; i++){
        std::vector<double> distances(n);
        std::vector<vertex_descriptor> p(n);
        vertex_descriptor from = boost::vertex(i, matching_graph);
        boost::dijkstra_shortest_paths(matching_graph, from,
            boost::weight_map(boost::get(&WeightedEdgeData::weight, matching_graph))
            .distance_map(boost::make_iterator_property_map(distances.begin(),
                            boost::get(boost::vertex_index, matching_graph)))
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


void MatchingGraph::ResetDijkstraNeighbours(){
    int n = boost::num_vertices(matching_graph);
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


std::vector<std::pair<int, double>> MatchingGraph::GetNearestNeighbours(
    int source, int num_neighbours, std::vector<int>& defect_id){
    int n = boost::num_vertices(matching_graph);
    if (source < 0 || source >= n) {
        throw std::invalid_argument("source must be non-negative and less than the number of nodes "
                                    "in the matching graph");
    }
    if (defect_id.size() != n) {
        throw std::invalid_argument("defect_id must have the same number of elements as the number "
                                    "of nodes in the matching graph");
    }
    if (num_neighbours < 0) {
        throw std::invalid_argument("num_neighbours must be a positive integer");
    }
    double inf = std::numeric_limits<double>::max();
    ResetDijkstraNeighbours();
    _distances[source] = 0;
    std::vector<int> examined_defects;
    std::vector<int> discovered_nodes;
    int source_is_defect = defect_id[source] > -1;
    DijkstraNeighbourVisitor vis = DijkstraNeighbourVisitor(
        defect_id, num_neighbours + source_is_defect, examined_defects, discovered_nodes);
    vertex_descriptor from = boost::vertex(source, matching_graph);
    try {
        boost::dijkstra_shortest_paths_no_color_map_no_init(matching_graph, from,
                &_predecessors[0], boost::make_iterator_property_map(_distances.begin(),
                boost::get(boost::vertex_index, matching_graph)),
                boost::get(&WeightedEdgeData::weight, matching_graph),
                boost::get(boost::vertex_index, matching_graph),
                std::less<double>(),
                boost::closed_plus<double>(),
                inf,
                0,
                vis);
    } catch (exit_search e) {}
    
    std::vector<std::pair<int, double>> neighbours;
    for (auto d : examined_defects){
        if (d != source){
            neighbours.push_back({d, _distances[d]});
        }
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


std::vector<int> MatchingGraph::GetPath(
    int source, int target){
    int n = boost::num_vertices(matching_graph);
    if (source >= n || target >= n
        || source < 0 || target < 0){
        throw std::invalid_argument("source and target must non-negative and less "
                                    "than the number of nodes");
    }
    double inf = std::numeric_limits<double>::max();
    ResetDijkstraNeighbours();
    _distances[source] = 0;
    std::vector<int> discovered_nodes;
    DijkstraPathVisitor vis = DijkstraPathVisitor(
        target, discovered_nodes);
    vertex_descriptor from = boost::vertex(source, matching_graph);
    try {
        boost::dijkstra_shortest_paths_no_color_map_no_init(matching_graph, from,
                &_predecessors[0], boost::make_iterator_property_map(_distances.begin(),
                boost::get(boost::vertex_index, matching_graph)),
                boost::get(&WeightedEdgeData::weight, matching_graph),
                boost::get(boost::vertex_index, matching_graph),
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


double MatchingGraph::Distance(int node1, int node2) {
    int num_nodes = GetNumNodes();
    if (node1 >= num_nodes || node2 >= num_nodes
        || node1 < 0 || node2 < 0){
        throw std::invalid_argument("node1 and node2 must non-negative and less "
                                    "than the number of nodes");
    }
    if (!HasComputedAllPairsShortestPaths()){
        ComputeAllPairsShortestPaths();
    }
    vertex_descriptor n2 = boost::vertex(node2, matching_graph);
    return all_distances[node1][n2];
}

std::vector<int> MatchingGraph::ShortestPath(int node1, int node2) {
    int num_nodes = GetNumNodes();
    if (node1 >= num_nodes || node2 >= num_nodes
        || node1 < 0 || node2 < 0){
        throw std::invalid_argument("node1 and node2 must non-negative and less "
                                    "than the number of nodes");
    }
    if (!HasComputedAllPairsShortestPaths()){
        ComputeAllPairsShortestPaths();
    }
    std::vector<vertex_descriptor> parent = all_predecessors[node2];
    auto index = boost::get(boost::vertex_index, matching_graph);
    int c = boost::vertex(node1, matching_graph);
    std::vector<int> path;
    path.push_back(index[c]);
    while (parent[c]!=c){
        c = parent[c];
        path.push_back(index[c]);
    }
    return path;
}


int MatchingGraph::GetNumEdges() const {
    return boost::num_edges(matching_graph);
}


int MatchingGraph::GetNumFaultIDs() const {
    auto qid = boost::get(&WeightedEdgeData::fault_ids, matching_graph);
    int num_edges = boost::num_edges(matching_graph);
    int maxid = -1;
    std::set<int> edge_fault_ids;
    auto es = boost::edges(matching_graph);
    for (auto eit = es.first; eit != es.second; ++eit) {
        edge_fault_ids = qid[*eit];
        for (auto fault_id : edge_fault_ids){
            if (fault_id >= maxid){
                maxid = fault_id;
            }
            if (fault_id < 0 && fault_id != -1){
                throw std::runtime_error("Fault ids must be non-negative, or -1 if no fault IDs are associated with the edge.");
            }
        }
    }
    return maxid + 1;
}

int MatchingGraph::GetNumNodes() const {
    return boost::num_vertices(matching_graph);
};

std::set<int> MatchingGraph::FaultIDs(int node1, int node2) const {
    int num_nodes = GetNumNodes();
    if (node1 >= num_nodes || node2 >= num_nodes
        || node1 < 0 || node2 < 0){
        throw std::invalid_argument("node1 and node2 must non-negative and less "
                                    "than the number of nodes");
    }
    auto e = boost::edge(node1, node2, matching_graph);
    if (!e.second){
        throw std::invalid_argument("Graph does not contain edge ("
                        + std::to_string((int)node1) + ", "
                        + std::to_string((int)node2) + ").");
    }
    return matching_graph[e.first].fault_ids;
}

std::pair<py::array_t<std::uint8_t>,py::array_t<std::uint8_t>> MatchingGraph::AddNoise() const {
    auto syndrome = new std::vector<int>(GetNumNodes(), 0);
    auto error = new std::vector<int>(GetNumFaultIDs(), 0);
    double p;
    std::set<int> qids;
    vertex_descriptor s, t;
    bool to_flip;
    auto es = boost::edges(matching_graph);
    for (auto eit = es.first; eit != es.second; ++eit) {
        p = matching_graph[*eit].error_probability;
        if ((p >= 0) && (rand_float(0.0, 1.0) < p)){
            s = boost::source(*eit, matching_graph);
            t = boost::target(*eit, matching_graph);
            (*syndrome)[s] = ((*syndrome)[s] + 1) % 2;
            (*syndrome)[t] = ((*syndrome)[t] + 1) % 2;
            qids = matching_graph[*eit].fault_ids;
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

std::set<int> MatchingGraph::GetBoundary() const {
    return boundary;
}

void MatchingGraph::SetBoundary(std::set<int>& boundary) {
    for (auto b: boundary){
        if (b < 0){
            throw std::invalid_argument("Boundary nodes must be non-negative.");
        }
    }
    this->boundary = boundary;
    connected_components_need_updating = true;
    return;
}

std::vector<std::tuple<int,int,WeightedEdgeData>> MatchingGraph::GetEdges() const {
    std::vector<std::tuple<int,int,WeightedEdgeData>> edges;
    auto es = boost::edges(matching_graph);
    for (auto eit = es.first; eit != es.second; ++eit) {
        WeightedEdgeData edata = matching_graph[*eit];
        int s = boost::source(*eit, matching_graph);
        int t = boost::target(*eit, matching_graph);
        if (edata.weight_is_negative) {
            edata.weight = -1 * edata.weight;
        }
        std::tuple<int,int,WeightedEdgeData> edge = std::make_tuple(s, t, edata);
        edges.push_back(edge);
    }
    return edges;
}

bool MatchingGraph::HasComputedAllPairsShortestPaths() const {
    int n = boost::num_vertices(matching_graph);
    bool has_distances = all_distances.size() == n;
    bool has_preds = all_predecessors.size() == n;
    return has_distances && has_preds;
}

int MatchingGraph::GetNumConnectedComponents() {
    if (connected_components_need_updating){
        component.resize(GetNumNodes());
        num_components = boost::connected_components(matching_graph, &component[0]);
        component_boundary.resize(num_components, -1);
        for (auto b : boundary){
            if (b >= GetNumNodes() || b < 0){
                throw std::invalid_argument(
                    "Boundary node ID " + std::to_string(b)
                    + " does not correspond to a node in the graph, which "
                    "has " + std::to_string(GetNumNodes()) + " nodes."
                );
            }
            int c = component[b];
            if (component_boundary[c] == -1){
                component_boundary[c] = b;
            }
        }
    }
    connected_components_need_updating = false;
    return num_components;
}

void MatchingGraph::FlipBoundaryNodesIfNeeded(std::set<int> &defects){
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

bool MatchingGraph::AllEdgesHaveErrorProbabilities() const {
    return all_edges_have_error_probabilities;
}

std::string MatchingGraph::repr() const {
    std::stringstream ss;
    ss << "<pymatching._cpp_mwpm.MatchingGraph object with ";
    ss << GetNumNodes() << " nodes, ";
    ss << GetNumEdges() << " edges and " << GetBoundary().size() << " boundary nodes>";
    return ss.str();
}
