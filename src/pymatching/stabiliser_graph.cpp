#include <cstdlib>
#include "stabiliser_graph.h"


double IStabiliserGraph::SpaceTimeDistance(int node1, int node2) const {
    int num_stab = GetNumStabilisers();
    if ((node1 < num_stab) && (node2 < num_stab)){
        return Distance(node1, node2);
    }
    int t1 = node1 / num_stab;
    int r1 = node1 % num_stab;
    int t2 = node2 / num_stab;
    int r2 = node2 % num_stab;
    return std::abs(t2-t1) + Distance(r1, r2);
}

std::vector<int> IStabiliserGraph::SpaceTimeShortestPath(int node1, int node2) const {
    int num_stab = GetNumStabilisers();
    int r1 = node1 % num_stab;
    int r2 = node2 % num_stab;
    return ShortestPath(r1, r2);
}