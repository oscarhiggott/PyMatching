#include <iostream>
#include "mwpm.h"
#include <map>
#include <algorithm>
#include <PerfectMatching.h>


void AddEdge(graph& g, int node1, int node2){
    g[node1].push_back(node2);
    g[node2].push_back(node1);
}


GraphData StabiliserGraph(const py::array_t<int>& indices, int num_stabilisers){
    auto x = indices.unchecked<1>();
    assert((x.shape(0) % 2)==0);
    graph g(num_stabilisers);
    EdgeData qubit;
    for (py::ssize_t i=0; i<(x.shape(0)/2); i++){
        AddEdge(g, x[2*i], x[2*i+1]);
        std::pair<int,int> edge = std::make_pair(std::min(x[2*i], x[2*i+1]),
                                        std::max(x[2*i], x[2*i+1]));
        qubit.insert({edge, i});
    }
    GraphData gd;
    gd.g = g;
    gd.qubit = qubit;
    return gd;
}


BFSResult BreadthFirstSearch(const graph& g, int source){
    std::vector<bool> visited(g.size(), false);
    std::list<int> queue;
    std::vector<int> distance(g.size(), INT_MAX);
    std::vector<int> parent(g.size(), -1);

    distance[source]=0;
    visited[source]=true;
    queue.push_back(source);

    while (!queue.empty()){
        int current = queue.front();
        queue.pop_front();
        for (vi_size_t i=0; i<g[current].size(); i++){
            int next = g[current][i];
            if (visited[next]==false){
                visited[next]=true;
                distance[next] = distance[current]+1;
                parent[next] = current;
                queue.push_back(next);
            }
        }
    }
    BFSResult res;
    res.distance = distance;
    res.parent = parent;
    return res;
}


APSPResult AllPairsShortestPath(const graph& g){
    vvint distances(g.size());
    vvint parents(g.size());
    for (graph::size_type i=0; i<g.size(); i++){
        auto res = BreadthFirstSearch(g, i);
        distances[i] = res.distance;
        parents[i] = res.parent;
    }
    APSPResult res;
    res.distances = distances;
    res.parents = parents;
    return res;
}


std::vector<int> ShortestPath(const std::vector<int>& parent, int dest){
    std::vector<int> path;
    int c = dest;
    path.push_back(c);
    while (parent[c]!=-1){
        c = parent[c];
        path.push_back(c);
    }
    return path;
}


py::array_t<int> Decode(const APSPResult& apsp, const py::array_t<int>& defects, 
                           const EdgeData& qubit, int num_qubits){
    auto d = defects.unchecked<1>();
    int num_nodes = d.shape(0);
    int num_edges = num_nodes*(num_nodes-1)/2;
    PerfectMatching *pm = new PerfectMatching(num_nodes, num_edges);
    for (py::size_t i = 0; i<num_nodes; i++){
        for (py::size_t j=i+1; j<num_nodes; j++){
            pm->AddEdge(i, j, apsp.distances[d(i)][d(j)]);
        }
    };
    pm->options.verbose = false;
    pm->Solve();
   auto correction = new std::vector<int>(num_qubits, 0);

    for (py::size_t i = 0; i<num_nodes; i++){
        int j = pm->GetMatch(i);
        if (i<j){
            std::vector<int> path = ShortestPath(apsp.parents[d(i)], d(j));
            for (vi_size_t k=0; k<path.size()-1; k++){
                int s1 = std::min(path[k],path[k+1]);
                int s2 = std::max(path[k],path[k+1]);
                int qid = qubit.find(std::make_pair(s1, s2))->second;
                (*correction)[qid] = 1;
            }
        }
    }
    delete pm;
    auto capsule = py::capsule(correction, [](void *correction) { delete reinterpret_cast<std::vector<int>*>(correction); });
    return py::array_t<int>(correction->size(), correction->data(), capsule);
}
