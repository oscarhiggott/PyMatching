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
        virtual int GetNumNodes() const;
        virtual int GetNumEdges() const;
        std::pair<py::array_t<std::uint8_t>,py::array_t<std::uint8_t>> AddNoise() const;
        std::vector<std::vector<double>> all_distances;
        std::vector<std::vector<vertex_descriptor>> all_predecessors;
        bool all_edges_have_error_probabilities;
        virtual std::vector<int> GetBoundary() const;
        virtual void SetBoundary(std::vector<int>& boundary);
        virtual bool HasComputedAllPairsShortestPaths() const;
        void ResetDijkstraNeighbours();
        std::vector<std::pair<int, double>> GetNearestNeighbours(
            int source, int num_neighbours, std::vector<int>& defect_id);
        std::vector<int> GetPath(
            int source, int target);
        virtual int GetNumConnectedComponents() const;
        std::vector<int> boundary;
        std::vector<double> _distances;
        std::vector<int> _predecessors;
};