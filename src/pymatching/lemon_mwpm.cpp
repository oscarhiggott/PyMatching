#include <iostream>
#include "lemon_mwpm.h"
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <vector>
#include <string>
#include "weighted_stabiliser_graph.h"
#include "stabiliser_graph.h"
#include <stdexcept>
#include <set>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <boost/graph/adjacency_list.hpp>

using namespace lemon;
using namespace std;

int LemonTest()
{
  ListDigraph g;
  ListDigraph::Node u = g.addNode();
  ListDigraph::Node v = g.addNode();
  ListDigraph::Arc  a = g.addArc(u, v);
  cout << "Hello World! This is LEMON library here." << endl;
  cout << "We have a directed graph with " << countNodes(g) << " nodes "
       << "and " << countArcs(g) << " arc." << endl;
  return 0;
}

double LemonMatchingTest()
{
    typedef lemon::ListGraph UGraph;
    typedef UGraph::EdgeMap<double> LengthMap;
    using lemon::INVALID;

    UGraph g;

    LengthMap length(g);

    UGraph::Node a = g.addNode();
    std::cout<<"a: "<<g.id(a)<<std::endl;
    UGraph::Node b = g.addNode();
    std::cout<<"b: "<<g.id(b)<<std::endl;
    UGraph::Node c = g.addNode();
    std::cout<<"c: "<<g.id(c)<<std::endl;
    UGraph::Node d = g.addNode();
    std::cout<<"d: "<<g.id(d)<<std::endl;
    UGraph::Edge e1 = g.addEdge(a,b);
    UGraph::Edge e2 = g.addEdge(b,c);
    UGraph::Edge e3 = g.addEdge(c,d);
    UGraph::Edge e4 = g.addEdge(d,a);

    length[e1] = 0.8;
    length[e2] = 0.9;
    length[e3] = 1.0;
    length[e4] = 1.1;

    typedef lemon::MaxWeightedPerfectMatching<UGraph,LengthMap> MWPM;
    MWPM gg(g, length);
    gg.run();
    std::cout<<"Matching: "<<gg.matchingWeight()<<std::endl;
    std::cout<<g.id(a)<<" with "<< g.id(gg.mate(a))<<std::endl;
    std::cout<<g.id(b)<<" with "<< g.id(gg.mate(b))<<std::endl;
    std::cout<<g.id(c)<<" with "<< g.id(gg.mate(c))<<std::endl;
    std::cout<<g.id(d)<<" with "<< g.id(gg.mate(d))<<std::endl;
    return gg.matchingWeight();
}

//typedef boost::property< boost::edge_weight_t, float, boost::property< boost::edge_index_t, int > >
//    Weight;
//typedef boost::adjacency_list< boost::vecS, boost::vecS, 
//        boost::undirectedS, boost::no_property, Weight >
//    defect_graph;

py::array_t<std::uint8_t> DecodeLemonMatchNeighbourhood(WeightedStabiliserGraph& sg, const py::array_t<int>& defects, int num_neighbours){
    auto d = defects.unchecked<1>();
    int num_defects = d.shape(0);

    int num_nodes = boost::num_vertices(sg.stabiliser_graph);

    std::vector<int> defect_id(num_nodes, -1);
    for (int i=0; i<num_defects; i++){
        defect_id[d(i)] = i;
    }

    num_neighbours = std::min(num_neighbours, num_defects-1) + 1;
    // int num_edges_max = (num_defects * num_neighbours);

    typedef lemon::ListGraph UGraph;
    typedef UGraph::EdgeMap<double> LengthMap;
    UGraph g;
    LengthMap length(g);
    UGraph::NodeMap<int> node_map(g);
    std::vector<UGraph::Node> node_list;
    for (int i=0; i<num_defects; i++){
        UGraph::Node x;
        x = g.addNode();
        node_map[x] = i;
        node_list.push_back(x);
    }

    // defect_graph g(num_defects);
    //std::vector< boost::graph_traits< defect_graph >::vertex_descriptor > match(num_defects);
    // PerfectMatching *pm = new PerfectMatching(num_defects, num_edges_max);
    std::vector<std::pair<int, double>> neighbours;
    std::vector<std::set<int>> adj_list(num_defects);
    int j;
    bool is_in;
    for (int i=0; i<num_defects; i++){
        neighbours = sg.GetNearestNeighbours(d(i), num_neighbours, defect_id);
        for (const auto &neighbour : neighbours){
            j = defect_id[neighbour.first];
            is_in = adj_list[i].find(j) != adj_list[i].end();
            if (!is_in && i!=j){
                UGraph::Edge e = g.addEdge(node_list[i], node_list[j]);
                length[e] = 1.0/(1.0+neighbour.second);
                // boost::add_edge(i, j, Weight(1.0/(1.0+neighbour.second)), g);
                adj_list[i].insert(j);
                adj_list[j].insert(i);
            }
        }
    }

//    boost::maximum_weighted_matching(g, &match[0]);
    typedef lemon::MaxWeightedPerfectMatching<UGraph,LengthMap> MWPM;
    MWPM pm(g, length);
    pm.run();

    int N = sg.GetNumQubits();
    auto correction = new std::vector<int>(N, 0);

    std::set<int> remaining_defects;
    for (int i=0; i<num_defects; i++){
        remaining_defects.insert(i);
    }

    std::vector<int> path;
    int i;
    std::set<int> qids;
    while (remaining_defects.size() > 0){
        i = *remaining_defects.begin();
        remaining_defects.erase(remaining_defects.begin());
        j = node_map[pm.mate(node_list[i])];
        remaining_defects.erase(j);
        path = sg.GetPath(d(i), d(j));
        for (std::vector<int>::size_type k=0; k<path.size()-1; k++){
            qids = sg.QubitIDs(path[k], path[k+1]);
            for (auto qid : qids){
                if ((qid != -1) && (qid >= 0) && (qid < N)){
                    (*correction)[qid] = ((*correction)[qid] + 1) % 2;
                }
            }
        }
    }
    auto capsule = py::capsule(correction, [](void *correction) { delete reinterpret_cast<std::vector<int>*>(correction); });
    return py::array_t<int>(correction->size(), correction->data(), capsule);
}