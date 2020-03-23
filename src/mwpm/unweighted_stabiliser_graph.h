#pragma once

#include <vector>
#include <list>
#include <climits>
#include <map>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "stabiliser_graph.h"

namespace py = pybind11;

typedef std::vector<std::vector<int>> graph;

typedef std::map<std::pair<int,int>, int> EdgeData;

typedef std::vector<int>::size_type vi_size_t;

struct BFSResult{
    std::vector<int> distance;
    std::vector<int> parent;
};

BFSResult BreadthFirstSearch(const graph& g, int source);

typedef std::vector<std::vector<int>> vvint;

struct APSPResult{
    vvint distances;
    vvint parents;
};

APSPResult AllPairsShortestPath(const graph& g);

std::vector<int> GetShortestPath(const std::vector<int>& parent, int dest);


class UnweightedStabiliserGraph : public IStabiliserGraph{
    public:
        graph adj_list;
        EdgeData qubit_ids;
        APSPResult shortest_paths;
        int num_qubits;
        int num_stabilisers;
        UnweightedStabiliserGraph(const py::array_t<int>& indices);
        void AddEdge(int node1, int node2, int qubit_id);
        virtual int Distance(int node1, int node2) const;
        // virtual int SpaceTimeDistance(int node1, int node2) const;
        virtual std::vector<int> ShortestPath(int node1, int node2) const;
        // virtual std::vector<int> SpaceTimeShortestPath(int node1, int node2) const;
        virtual int QubitID(int node1, int node2) const;
        virtual int GetNumQubits() const;
        virtual int GetNumStabilisers() const;
};