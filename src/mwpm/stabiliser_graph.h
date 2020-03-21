#pragma once

#include <vector>

class IStabiliserGraph{
    public:
        virtual int Distance(int node1, int node2) const = 0;
        virtual int SpaceTimeDistance(int node1, int node2) const = 0;
        virtual std::vector<int> ShortestPath(int node1, int node2) const = 0;
        virtual std::vector<int> SpaceTimeShortestPath(int node1, int node2) const = 0;
        virtual int QubitID(int node1, int node2) const = 0;
        virtual int GetNumQubits() const = 0;
};