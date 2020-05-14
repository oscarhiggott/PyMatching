#include <map>
#include <algorithm>
#include "unweighted_stabiliser_graph.h"


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


std::vector<int> GetShortestPath(const std::vector<int>& parent, int dest){
    std::vector<int> path;
    int c = dest;
    path.push_back(c);
    while (parent[c]!=-1){
        c = parent[c];
        path.push_back(c);
    }
    return path;
}

UnweightedStabiliserGraph::UnweightedStabiliserGraph(
    int num_stabilisers, int boundary
) : num_stabilisers(num_stabilisers), num_qubits(0), boundary(boundary) {
    adj_list.resize(num_stabilisers);
}

UnweightedStabiliserGraph::UnweightedStabiliserGraph(
    const py::array_t<int>& indices, int boundary
    ) : boundary(boundary){
        auto x = indices.unchecked<1>();
        assert((x.shape(0) % 2)==0);
        int smax = 0;
        for (py::ssize_t i=0; i<(x.shape(0)); i++){
            if (x[i] > smax){
                smax = x[i];
            }
        }
        this->num_stabilisers = smax+1;
        this->num_qubits = 0;
        adj_list.resize(num_stabilisers);
        for (py::ssize_t i=0; i<(x.shape(0)/2); i++){
            AddEdge(x[2*i], x[2*i+1], i);
        }
    shortest_paths = AllPairsShortestPath(adj_list);
}


void UnweightedStabiliserGraph::AddEdge(int node1, int node2, int qubit_id){
    assert(node1 < num_stabilisers);
    assert(node2 < num_stabilisers);
    adj_list[node1].push_back(node2);
    adj_list[node2].push_back(node1);
    std::pair<int,int> edge = std::make_pair(std::min(node1, node2),
                                    std::max(node1, node2));
    qubit_ids.insert({edge, qubit_id});
    num_qubits++;
}


double UnweightedStabiliserGraph::Distance(int node1, int node2) const {
    assert(node1 < num_stabilisers);
    assert(node2 < num_stabilisers);
    return (double) shortest_paths.distances[node1][node2];
}

std::vector<int> UnweightedStabiliserGraph::ShortestPath(int node1, int node2) const {
    assert(node1 < num_stabilisers);
    assert(node2 < num_stabilisers);
    return GetShortestPath(shortest_paths.parents[node2], node1);
}

int UnweightedStabiliserGraph::QubitID(int node1, int node2) const {
    int s1 = std::min(node1, node2);
    int s2 = std::max(node1, node2);
    int qid = qubit_ids.find(std::make_pair(s1, s2))->second;
    assert(qid < num_qubits);
    return qid;
}

int UnweightedStabiliserGraph::GetNumQubits() const {
    return num_qubits;
}

int UnweightedStabiliserGraph::GetNumStabilisers() const {
    return num_stabilisers;
}

void UnweightedStabiliserGraph::ComputeAllPairsShortestPaths(){
    shortest_paths = AllPairsShortestPath(adj_list);
}

int UnweightedStabiliserGraph::GetBoundaryVertex() const {
    return boundary;
}

void UnweightedStabiliserGraph::SetBoundaryVertex(int boundary) {
    this->boundary = boundary;
    return;
}