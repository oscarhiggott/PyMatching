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
        virtual int GetNumQubits() const = 0;
        virtual int GetNumStabilisers() const = 0;
        virtual std::vector<int> GetBoundary() const = 0;
        virtual void SetBoundary(std::vector<int>& boundary) = 0;
        virtual bool HasComputedAllPairsShortestPaths() const = 0;
        virtual void ComputeAllPairsShortestPaths() = 0;
};