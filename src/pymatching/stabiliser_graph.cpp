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

#include <cstdlib>
#include "stabiliser_graph.h"


double IStabiliserGraph::SpaceTimeDistance(int node1, int node2) {
    int num_stab = GetNumNodes();
    if ((node1 < num_stab) && (node2 < num_stab)){
        return Distance(node1, node2);
    }
    int t1 = node1 / num_stab;
    int r1 = node1 % num_stab;
    int t2 = node2 / num_stab;
    int r2 = node2 % num_stab;
    return std::abs(t2-t1) + Distance(r1, r2);
}

std::vector<int> IStabiliserGraph::SpaceTimeShortestPath(int node1, int node2) {
    int num_stab = GetNumNodes();
    int r1 = node1 % num_stab;
    int r2 = node2 % num_stab;
    return ShortestPath(r1, r2);
}