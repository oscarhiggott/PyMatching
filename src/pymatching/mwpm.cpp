#include "mwpm.h"
#include <PerfectMatching.h>
#include "stabiliser_graph.h"
#include "weighted_stabiliser_graph.h"
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <ctime>
#include <set>


py::array_t<std::uint8_t> Decode(IStabiliserGraph& sg, const py::array_t<int>& defects){
    if (!sg.HasComputedAllPairsShortestPaths()){
        sg.ComputeAllPairsShortestPaths();
    }
    auto d = defects.unchecked<1>();
    int num_nodes = d.shape(0);
    int num_edges = num_nodes*(num_nodes-1)/2;

    PerfectMatching *pm = new PerfectMatching(num_nodes, num_edges);
    for (py::size_t i = 0; i<num_nodes; i++){
        for (py::size_t j=i+1; j<num_nodes; j++){
            pm->AddEdge(i, j, sg.SpaceTimeDistance(d(i), d(j)));
        }
    };
    pm->options.verbose = false;
    pm->Solve();
    int N = sg.GetNumQubits();
    auto correction = new std::vector<int>(N, 0);
    std::set<int> qids;
    for (py::size_t i = 0; i<num_nodes; i++){
        int j = pm->GetMatch(i);
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
    delete pm;

    auto capsule = py::capsule(correction, [](void *correction) { delete reinterpret_cast<std::vector<int>*>(correction); });
    return py::array_t<int>(correction->size(), correction->data(), capsule);
}


py::array_t<std::uint8_t> DecodeMatchNeighbourhood(WeightedStabiliserGraph& sg, const py::array_t<int>& defects, int num_neighbours){
    auto d = defects.unchecked<1>();
    int num_defects = d.shape(0);

    int num_nodes = boost::num_vertices(sg.stabiliser_graph);

    std::vector<int> defect_id(num_nodes, -1);
    for (int i=0; i<num_defects; i++){
        defect_id[d(i)] = i;
    }

    num_neighbours = std::min(num_neighbours, num_defects-1) + 1;
    int num_edges_max = (num_defects * num_neighbours);
    PerfectMatching *pm = new PerfectMatching(num_defects, num_edges_max);
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
                pm->AddEdge(i, j, neighbour.second);
                adj_list[i].insert(j);
                adj_list[j].insert(i);
            }
        }
    }

    pm->options.verbose = false;
    pm->Solve();

    int N = sg.GetNumQubits();
    auto correction = new std::vector<int>(N, 0);

    std::set<int> remaining_defects;
    for (int i=0; i<num_defects; i++){
        remaining_defects.insert(i);
    }

    std::vector<int> path;
    int i;
    std::set<int> qids;
    while (remaining_defects.size() > 0){
        i = *remaining_defects.begin();
        remaining_defects.erase(remaining_defects.begin());
        j = pm->GetMatch(i);
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
    delete pm;
    auto capsule = py::capsule(correction, [](void *correction) { delete reinterpret_cast<std::vector<int>*>(correction); });
    return py::array_t<int>(correction->size(), correction->data(), capsule);
}