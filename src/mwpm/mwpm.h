#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <list>
#include <climits>
#include <map>

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

class StabiliserGraph{
    public:
    graph adj_list;
    EdgeData qubit_ids;
    APSPResult shortest_paths;
    int num_qubits;
    int num_stabilisers;
    StabiliserGraph(const py::array_t<int>& indices, int num_stabilisers, int num_qubits);
    void AddEdge(int node1, int node2);
    int Distance(int node1, int node2) const;
    int SpaceTimeDistance(int node1, int node2) const;
    std::vector<int> ShortestPath(int node1, int node2) const;
    std::vector<int> SpaceTimeShortestPath(int node1, int node2) const;
    int QubitID(int node1, int node2) const;
};

py::array_t<int> Decode(const StabiliserGraph& sg, const py::array_t<int>& defects);
