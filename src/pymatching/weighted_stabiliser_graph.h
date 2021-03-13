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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <memory>
#include <set>
#include <cstdint>
#include "stabiliser_graph.h"

namespace py = pybind11;

struct WeightedEdgeData {
    std::set<int> qubit_ids;
    double weight;
    double error_probability;
    bool has_error_probability;
};

typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
    boost::no_property, WeightedEdgeData > wgraph_t;
typedef boost::graph_traits < wgraph_t >::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits < wgraph_t >::edge_descriptor edge_descriptor;


/**
 * @brief 
 * 
 */
class WeightedStabiliserGraph : public IStabiliserGraph{
    public:
        wgraph_t stabiliser_graph;
        /**
         * @brief Construct a new Weighted Stabiliser Graph object
         * 
         * @param num_stabilisers Number of stabiliser nodes (this excludes boundary nodes)
         * @param boundary Indices of the boundary nodes
         */
        WeightedStabiliserGraph(
            int num_stabilisers,
            std::vector<int>& boundary
            );
        /**
         * @brief Add an edge to the Weighted Stabiliser Graph object
         * 
         * @param node1 Index of the first node
         * @param node2 Index of the second node
         * @param qubit_ids Indices of the qubits associated with this edge
         * @param weight Weight of the edge
         * @param error_probability The error probability associated with an edge (optional, set to -1 if not needed)
         * @param has_error_probability Flag whether a valid error probability has been supplied
         */
        void AddEdge(int node1, int node2, std::set<int> qubit_ids,
                     double weight, double error_probability,
                     bool has_error_probability);
        /**
         * @brief Compute and store the shortest path between every pair of nodes in the matching graph.
         * This is only used or needed if exact matching is used (rather than the default local matching).
         * Note that this is only useful for relatively small matching graphs, and is very memory and 
         * compute intensive matching graphs with many thousands of nodes.
         * 
         */
        virtual void ComputeAllPairsShortestPaths();
        /**
         * @brief Distance between two nodes in the matching graph. This method is used only for exact matching,
         * since it uses the pre-computed all-pairs shortest paths data, and computes the all-pairs shortest paths
         * if this has not yet been done.
         * 
         * @param node1 The index of the first node
         * @param node2 The index of the second node
         * @return double The distance between the first node and the second node in the matching graph
         */
        virtual double Distance(int node1, int node2);
        /**
         * @brief Shortest path between two nodes in the matching graph. This method is used only for exact matching,
         * since it uses the pre-computed all-pairs shortest paths data, and computes the all-pairs shortest paths
         * if this has not yet been done.
         * 
         * @param node1 The index of the first node
         * @param node2 The index of the second node
         * @return std::vector<int> The vertex indices along the path between the first node and the second node in the matching graph
         */
        virtual std::vector<int> ShortestPath(int node1, int node2);
        /**
         * @brief Get the qubit ids associated with the edge between node1 and node2
         * 
         * @param node1 The index of the first node
         * @param node2 The index of the second node
         * @return std::set<int> The set of qubit ids associated with the edge between node1 and node2
         */
        virtual std::set<int> QubitIDs(int node1, int node2) const;
        /**
         * @brief Get the number of qubits associated with edges of the matching graph
         * 
         * @return int The number of qubits
         */
        virtual int GetNumQubits() const;
        /**
         * @brief Get the number of nodes in the stabiliser graph (this includes both boundaries and stabilisers)
         * 
         * @return int Number of nodes in stabiliser_graph
         */
        virtual int GetNumNodes() const;
        /**
         * @brief Get the number of edges in the stabiliser graph
         * 
         * @return int Number of edges in stabiliser_graph
         */
        virtual int GetNumEdges() const;
        /**
         * @brief If an error_probability is assigned to every edge, flip each edge 
         * with its corresponding error_probability. If an edge is flipped, add (mod 2) 
         * an error to the associated qubits (specified by the qubit_ids edge data).
         * The qubit errors are returned as a binary numpy array, as well as a syndrome 
         * vector, also as a binary numpy array. The length of the syndrome vector is 
         * is the number of nodes in the stabiliser graph (there is an element for each 
         * stabiliser as well as for each boundary node). The syndromes of the boundary 
         * nodes are all set to zero unless the parity of the stabiliser syndromes is odd, 
         * in which case the first boundary node's syndrome is flipped.
         * 
         * @return std::pair<py::array_t<std::uint8_t>,py::array_t<std::uint8_t>> The first item in 
         * the pair is the noise vector, and the second item is the syndrome vector.
         */
        std::pair<py::array_t<std::uint8_t>,py::array_t<std::uint8_t>> AddNoise() const;
        /**
         * @brief The distance between every pair of nodes in the stabiliser graph, if 
         * ComputeAllPairsShortestPaths has been run. all_distances[i][j] is the distance 
         * between node i and node j in the stabiliser graph. Note that this is only used 
         * if exact matching is used (the function LemonDecode in lemon_mwpm.cpp), and not if 
         * (the Python default) local matching is used (the function LemonDecodeMatchNeighbourhood 
         * in lemon_mwpm.cpp).
         * 
         */
        std::vector<std::vector<double>> all_distances;
        /**
         * @brief all_predecessors[i] is the PredecessorMap found by Boost's dijkstra_shortest_paths 
         * from node i in stabiliser_graph to all other nodes. In other words, all_predecessors[i][j] 
         * is the parent of node j in the shortest path from node i to node j in stabiliser_graph.
         * Note that this is only used if exact matching is used (the function LemonDecode in lemon_mwpm.cpp), 
         * and not if (the Python default) local matching is used (the function LemonDecodeMatchNeighbourhood 
         * in lemon_mwpm.cpp).
         * 
         */
        std::vector<std::vector<vertex_descriptor>> all_predecessors;
        bool all_edges_have_error_probabilities; // true if every edge in stabiliser_graph has an error_probability in its edge data
        /**
         * @brief Get the indices of the boundary nodes
         * 
         * @return std::vector<int> The indices of the boundary nodes
         */
        virtual std::vector<int> GetBoundary() const;
        /**
         * @brief Get the edges of the stabiliser graph and their edge data
         * 
         * @return std::vector<std::tuple<int,int,WeightedEdgeData>> 
         */
        std::vector<std::tuple<int,int,WeightedEdgeData>> GetEdges() const;
        /**
         * @brief Set the indices of the boundary nodes
         * 
         * @param boundary The indices of the boundary nodes
         */
        virtual void SetBoundary(std::vector<int>& boundary);
        /**
         * @brief Flag whether or not the all-pairs shortest paths have been computed. 
         * The all-pairs shortest paths are only needed for exact matching, which is not 
         * a default option in the Python bindings.
         * 
         * @return true 
         * @return false 
         */
        virtual bool HasComputedAllPairsShortestPaths() const;
        /**
         * @brief Reset the _predecessors and _distances attributes, which 
         * are used by WeightedStabiliserGraph::GetPath and WeightedStabiliserGraph::GetNearestNeighbours. 
         * This just preallocates both of these vectors appropriately for use by the local dijkstra search 
         * and will only do so if these vectors are not already the correct size. Once these attributes 
         * are correctly preallocated the first time, the GetPath and GetNearestNeighbours methods always reset 
         * only the elements in these vectors that have been used for efficiency (since these methods only use 
         * a local search).
         * 
         */
        void ResetDijkstraNeighbours();
        /**
         * @brief Get the num_neighbours nearest neighbours i of a source node in the stabiliser_graph for which 
         * defect_id[i] > -1. This is a local Dijkstra search that halts once num_neighbours nodes i have been found 
         * that satisfy defect_id[i] > -1. The function returns a vector of pairs, where the first item in each 
         * pair is the distance from source to one of the nearest nodes i, and the second item in the pair 
         * is defect_id[i].
         * This is used by LemonDecodeMatchNeighbourhood to find the num_neigbours nearest defects to be included 
         * in the matching graph.
         * 
         * @param source 
         * @param num_neighbours 
         * @param defect_id 
         * @return std::vector<std::pair<int, double>> 
         */
        std::vector<std::pair<int, double>> GetNearestNeighbours(
            int source, int num_neighbours, std::vector<int>& defect_id);
        /**
         * @brief Get the shortest path from source to target using a version of Dijkstra 
         * that terminates once target is found.
         * 
         * @param source 
         * @param target 
         * @return std::vector<int> 
         */
        std::vector<int> GetPath(
            int source, int target);
        /**
         * @brief Get the number of connected components in the stabiliser graph
         * 
         * @return int The number of components
         */
        virtual int GetNumConnectedComponents() const;
        /**
         * @brief The indices of the boundary nodes
         * 
         */
        std::vector<int> boundary;
        std::vector<double> _distances;
        std::vector<int> _predecessors;
};