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
#include <pybind11/stl.h>
#include <memory>
#include <set>
#include <string>
#include <sstream>
#include <cstdint>

namespace py = pybind11;

struct WeightedEdgeData {
    std::set<int> fault_ids;
    double weight;
    double error_probability;
    bool has_error_probability;
    bool weight_is_negative;
    WeightedEdgeData();
    WeightedEdgeData(
        std::set<int> fault_ids,
        double weight,
        double error_probability,
        bool has_error_probability,
        bool weight_is_negative
    );
    std::string repr() const;
};

typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
    boost::no_property, WeightedEdgeData > wgraph_t;
typedef boost::graph_traits < wgraph_t >::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits < wgraph_t >::edge_descriptor edge_descriptor;


/**
 * @brief 
 * 
 */
class MatchingGraph{
    public:
        /**
         * @brief Construct a new Matching Graph object
         * 
         */
        MatchingGraph();
        /**
         * @brief Construct a new Matching Graph object
         * 
         * @param num_detectors Number of detector nodes (this excludes boundary nodes)
         * @param boundary Indices of the boundary nodes
         */
        MatchingGraph(
            int num_detectors,
            std::set<int>& boundary
            );
        /**
         * @brief Add an edge to the Matching Graph object
         * 
         * @param node1 Index of the first node
         * @param node2 Index of the second node
         * @param fault_ids Indices of the faults associated with this edge
         * @param weight Weight of the edge
         * @param error_probability The error probability associated with an edge (optional, set to -1 if not needed)
         * @param has_error_probability Flag whether a valid error probability has been supplied
         */
        void AddEdge(int node1, int node2, std::set<int> fault_ids,
                     double weight, double error_probability,
                     bool has_error_probability);
        /**
         * @brief Compute and store the shortest path between every pair of nodes in the matching graph.
         * This is only used or needed if exact matching is used (rather than the default local matching).
         * Note that this is only useful for relatively small matching graphs, and is very memory and 
         * compute intensive matching graphs with many thousands of nodes.
         * 
         */
        void ComputeAllPairsShortestPaths();
        /**
         * @brief Distance between two nodes in the matching graph. This method is used only for exact matching,
         * since it uses the pre-computed all-pairs shortest paths data, and computes the all-pairs shortest paths
         * if this has not yet been done.
         * 
         * @param node1 The index of the first node
         * @param node2 The index of the second node
         * @return double The distance between the first node and the second node in the matching graph
         */
        double Distance(int node1, int node2);
        /**
         * @brief Shortest path between two nodes in the matching graph. This method is used only for exact matching,
         * since it uses the pre-computed all-pairs shortest paths data, and computes the all-pairs shortest paths
         * if this has not yet been done.
         * 
         * @param node1 The index of the first node
         * @param node2 The index of the second node
         * @return std::vector<int> The vertex indices along the path between the first node and the second node in the matching graph
         */
        std::vector<int> ShortestPath(int node1, int node2);
        /**
         * @brief Get the fault ids associated with the edge between node1 and node2
         * 
         * @param node1 The index of the first node
         * @param node2 The index of the second node
         * @return std::set<int> The set of fault ids associated with the edge between node1 and node2
         */
        std::set<int> FaultIDs(int node1, int node2) const;
        /**
         * @brief Get the number of fault ids associated with edges of the matching graph
         * 
         * @return int The number of fault IDs
         */
        int GetNumFaultIDs() const;
        /**
         * @brief Get the number of nodes in the matching graph (this includes both boundaries and stabilisers)
         * 
         * @return int Number of nodes in matching_graph
         */
        int GetNumNodes() const;
        /**
         * @brief Get the number of edges in the matching graph
         * 
         * @return int Number of edges in matching_graph
         */
        int GetNumEdges() const;
        /**
         * @brief If an error_probability is assigned to every edge, flip each edge 
         * with its corresponding error_probability. If an edge is flipped, flip the corresponding
         * self-inverse faults (specified by the fault_ids edge data).
         * The fault ids to be flipped for the recovery are returned as a binary numpy array, as well as a syndrome
         * vector, also as a binary numpy array. The length of the syndrome vector is 
         * is the number of nodes in the matching graph (there is an element for each
         * stabiliser as well as for each boundary node). The syndromes of the boundary 
         * nodes are all set to zero unless the parity of the stabiliser syndromes is odd, 
         * in which case the first boundary node's syndrome is flipped.
         * 
         * @return std::pair<py::array_t<std::uint8_t>,py::array_t<std::uint8_t>> The first item in 
         * the pair is the noise vector, and the second item is the syndrome vector.
         */
        std::pair<py::array_t<std::uint8_t>,py::array_t<std::uint8_t>> AddNoise() const;
        /**
         * @brief Get the indices of the boundary nodes
         * 
         * @return std::vector<int> The indices of the boundary nodes
         */
        std::set<int> GetBoundary() const;
        /**
         * @brief Get the edges of the matching graph and their edge data
         * 
         * @return std::vector<std::tuple<int,int,WeightedEdgeData>> 
         */
        std::vector<std::tuple<int,int,WeightedEdgeData>> GetEdges() const;
        /**
         * @brief Set the indices of the boundary nodes
         * 
         * @param boundary The indices of the boundary nodes
         */
        void SetBoundary(std::set<int>& boundary);
        /**
         * @brief Flag whether or not the all-pairs shortest paths have been computed. 
         * The all-pairs shortest paths are only needed for exact matching, which is not 
         * a default option in the Python bindings.
         * 
         * @return true 
         * @return false 
         */
        bool HasComputedAllPairsShortestPaths() const;
        /**
         * @brief Reset the _predecessors and _distances attributes, which 
         * are used by MatchingGraph::GetPath and MatchingGraph::GetNearestNeighbours.
         * This just preallocates both of these vectors appropriately for use by the local dijkstra search 
         * and will only do so if these vectors are not already the correct size. Once these attributes 
         * are correctly preallocated the first time, the GetPath and GetNearestNeighbours methods always reset 
         * only the elements in these vectors that have been used for efficiency (since these methods only use 
         * a local search).
         * 
         */
        void ResetDijkstraNeighbours();
        /**
         * @brief Get the num_neighbours nearest neighbours i of a source node in the matching_graph for which
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
         * @brief Get the number of connected components in the matching graph
         * 
         * @return int The number of components
         */
        int GetNumConnectedComponents();
        /**
        * @brief Whether or not all edges in the graph are assigned error probabilities
        *
        * @ return bool True if all edges have error probabilities, else False
        */
        bool AllEdgesHaveErrorProbabilities() const;
        void FlipBoundaryNodesIfNeeded(std::set<int> &defects);
        std::string repr() const;
        void HandleNewNegativeWeightEdge(int u, int v, double weight, std::set<int> &fault_ids);
        std::set<int> negative_edge_syndrome;
        std::set<int> negative_edge_fault_ids;
        double negative_weight_sum;
    private:
        /**
         * @brief The indices of the boundary nodes
         *
         */
        std::set<int> boundary;
        std::vector<double> _distances;
        std::vector<int> _predecessors;
        int num_components;
        std::vector<int> component;
        std::vector<int> component_boundary;
        bool connected_components_need_updating;
        // true if every edge in matching_graph has an error_probability in its edge data
        bool all_edges_have_error_probabilities;
        wgraph_t matching_graph;
        /**
         * @brief The distance between every pair of nodes in the matching graph, if
         * ComputeAllPairsShortestPaths has been run. all_distances[i][j] is the distance
         * between node i and node j in the matching graph. Note that this is only used
         * if exact matching is used (the function LemonDecode in lemon_mwpm.cpp), and not if
         * (the Python default) local matching is used (the function LemonDecodeMatchNeighbourhood
         * in lemon_mwpm.cpp).
         *
         */
        std::vector<std::vector<double>> all_distances;
        /**
         * @brief all_predecessors[i] is the PredecessorMap found by Boost's dijkstra_shortest_paths
         * from node i in matching_graph to all other nodes. In other words, all_predecessors[i][j]
         * is the parent of node j in the shortest path from node i to node j in matching_graph.
         * Note that this is only used if exact matching is used (the function LemonDecode in lemon_mwpm.cpp),
         * and not if (the Python default) local matching is used (the function LemonDecodeMatchNeighbourhood
         * in lemon_mwpm.cpp).
         *
         */
        std::vector<std::vector<vertex_descriptor>> all_predecessors;
};