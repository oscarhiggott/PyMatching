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

#include <iostream>
#include "lemon_mwpm.h"
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <lemon/connectivity.h>
#include <vector>
#include <string>
#include "weighted_stabiliser_graph.h"
#include "stabiliser_graph.h"
#include <stdexcept>
#include <set>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <boost/graph/adjacency_list.hpp>

typedef lemon::ListGraph UGraph;
typedef UGraph::EdgeMap<double> LengthMap;


class DefectGraph {
    public:
        DefectGraph(int num_nodes);
        void AddEdge(int i, int j, double weight);
        UGraph g;
        LengthMap length;
        UGraph::NodeMap<int> node_map;
        std::vector<UGraph::Node> node_list;
        int num_nodes;
};

DefectGraph::DefectGraph(int num_nodes) : num_nodes(num_nodes),
         length(g), node_map(g)
{
    for (int i=0; i<num_nodes; i++){
        UGraph::Node x;
        x = g.addNode();
        node_map[x] = i;
        node_list.push_back(x);
    }
}

void DefectGraph::AddEdge(int i, int j, double weight){
    UGraph::Edge e = g.addEdge(node_list[i], node_list[j]);
    length[e] = weight;
}


MatchingResult LemonDecode(IStabiliserGraph& sg, const py::array_t<int>& defects, bool return_weight){
    MatchingResult matching_result;
    if (!sg.HasComputedAllPairsShortestPaths()){
        sg.ComputeAllPairsShortestPaths();
    }
    auto d = defects.unchecked<1>();
    int num_nodes = d.shape(0);

    DefectGraph defect_graph(num_nodes);

    for (py::size_t i = 0; i<num_nodes; i++){
        for (py::size_t j=i+1; j<num_nodes; j++){
            defect_graph.AddEdge(i, j, -1.0*sg.SpaceTimeDistance(d(i), d(j)));
        }
    };
    typedef lemon::MaxWeightedPerfectMatching<UGraph,LengthMap> MWPM;
    MWPM pm(defect_graph.g, defect_graph.length);
    pm.run();

    if (return_weight) {
        matching_result.weight = -1*pm.matchingWeight();
    } else {
        matching_result.weight = -1.0;
    }

    int N = sg.GetNumQubits();
    auto correction = new std::vector<int>(N, 0);
    std::set<int> qids;
    for (py::size_t i = 0; i<num_nodes; i++){
        int j = defect_graph.node_map[pm.mate(defect_graph.node_list[i])];
        if (i<j){
            std::vector<int> path = sg.SpaceTimeShortestPath(d(i), d(j));
            for (std::vector<int>::size_type k=0; k<path.size()-1; k++){
                qids = sg.QubitIDs(path[k], path[k+1]);
                for (auto qid : qids){
                    if ((qid != -1) && (qid >= 0) && (qid < N)){
                        (*correction)[qid] = ((*correction)[qid] + 1) % 2;
                    }
                }
            }
        }
    }

    auto capsule = py::capsule(correction, [](void *correction) { delete reinterpret_cast<std::vector<int>*>(correction); });
    auto corr = py::array_t<int>(correction->size(), correction->data(), capsule);
    matching_result.correction = corr;
    return matching_result;
}


MatchingResult LemonDecodeMatchNeighbourhood(WeightedStabiliserGraph& sg, const py::array_t<int>& defects, 
                                             int num_neighbours, bool return_weight){
    MatchingResult matching_result;
    auto d = defects.unchecked<1>();
    int num_defects = d.shape(0);

    int num_nodes = boost::num_vertices(sg.stabiliser_graph);

    std::vector<int> defect_id(num_nodes, -1);
    for (int i=0; i<num_defects; i++){
        defect_id[d(i)] = i;
    }
    num_neighbours = std::min(num_neighbours, num_defects-1) + 1;
    --num_neighbours;
    bool is_connected = false;
    std::unique_ptr<DefectGraph> defect_graph;
    while (!is_connected && num_neighbours < num_defects){
        ++num_neighbours;
        defect_graph = std::make_unique<DefectGraph>(num_defects);
        std::vector<std::pair<int, double>> neighbours;
        std::vector<std::set<int>> adj_list(num_defects);
        int j;
        bool is_in;
        for (int i=0; i<num_defects; i++){
            neighbours = sg.GetNearestNeighbours(d(i), num_neighbours, defect_id);
            for (const auto &neighbour : neighbours){
                j = defect_id[neighbour.first];
                is_in = adj_list[i].find(j) != adj_list[i].end();
                if (!is_in && i!=j){
                    defect_graph->AddEdge(i, j, -1.0*neighbour.second);
                    adj_list[i].insert(j);
                    adj_list[j].insert(i);
                }
            }
        }
        is_connected = lemon::connected(defect_graph->g);
    }

    if (!is_connected){
        throw std::runtime_error("Graph must have only one connected component");
    }

    typedef lemon::MaxWeightedPerfectMatching<UGraph,LengthMap> MWPM;
    MWPM pm(defect_graph->g, defect_graph->length);
    pm.run();

    if (return_weight) {
        matching_result.weight = -1*pm.matchingWeight();
    } else {
        matching_result.weight = -1.0;
    }

    int N = sg.GetNumQubits();
    auto correction = new std::vector<int>(N, 0);

    std::set<int> remaining_defects;
    for (int i=0; i<num_defects; i++){
        remaining_defects.insert(i);
    }

    std::vector<int> path;
    int i;
    int j;
    std::set<int> qids;
    while (remaining_defects.size() > 0){
        i = *remaining_defects.begin();
        remaining_defects.erase(remaining_defects.begin());
        j = defect_graph->node_map[pm.mate(defect_graph->node_list[i])];
        remaining_defects.erase(j);
        path = sg.GetPath(d(i), d(j));
        for (std::vector<int>::size_type k=0; k<path.size()-1; k++){
            qids = sg.QubitIDs(path[k], path[k+1]);
            for (auto qid : qids){
                if ((qid != -1) && (qid >= 0) && (qid < N)){
                    (*correction)[qid] = ((*correction)[qid] + 1) % 2;
                }
            }
        }
    }
    auto capsule = py::capsule(correction, [](void *correction) { delete reinterpret_cast<std::vector<int>*>(correction); });
    auto corr = py::array_t<int>(correction->size(), correction->data(), capsule);
    
    matching_result.correction = corr;
    return matching_result;
}