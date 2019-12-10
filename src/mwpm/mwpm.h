#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <list>
#include <climits>
#include <map>

namespace py = pybind11;

typedef std::vector<std::vector<int>> graph;

void AddEdge(graph& g, int node1, int node2);

typedef std::map<std::pair<int,int>, int> EdgeData;

struct GraphData{
    graph g;
    EdgeData qubit;
};

GraphData StabiliserGraph(const py::array_t<int>& indices, int num_stabilisers);

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

class StabiliserGraphCls{
    public:
    graph adj_list;
    EdgeData qubit_ids;
    APSPResult shortest_paths;
    int num_qubits;
    int num_stabilisers;
    StabiliserGraphCls(const py::array_t<int>& indices, int num_stabilisers, int num_qubits);
    int Distance(int i, int j) const;
    std::vector<int> ShortestPath(int i, int j) const;
    int QubitID(int i, int j) const;
};

py::array_t<int> Decode(const StabiliserGraphCls& sg, const py::array_t<int>& defects);
