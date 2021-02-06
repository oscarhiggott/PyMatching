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

#pragma once

#include <vector>
#include <set>

class IStabiliserGraph{
    public:
        virtual double Distance(int node1, int node2) = 0;
        virtual double SpaceTimeDistance(int node1, int node2);
        virtual std::vector<int> ShortestPath(int node1, int node2) = 0;
        virtual std::vector<int> SpaceTimeShortestPath(int node1, int node2);
        virtual std::set<int> QubitIDs(int node1, int node2) const = 0;
        virtual int GetNumEdges() const = 0;
        virtual int GetNumQubits() const = 0;
        virtual int GetNumNodes() const = 0;
        virtual int GetNumConnectedComponents() const = 0;
        virtual std::vector<int> GetBoundary() const = 0;
        virtual void SetBoundary(std::vector<int>& boundary) = 0;
        virtual bool HasComputedAllPairsShortestPaths() const = 0;
        virtual void ComputeAllPairsShortestPaths() = 0;
};