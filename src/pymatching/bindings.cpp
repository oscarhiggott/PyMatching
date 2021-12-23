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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "matching_graph.h"
#include "rand_gen.h"
#include "lemon_mwpm.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(_cpp_mwpm, m) {
     py::class_<WeightedEdgeData>(m, "WeightedEdgeData")
     .def(py::init<>(), u8R"(
     Initialises a WeightedEdgeData object
     )")
     .def(py::init<std::set<int>, double, double, bool, bool>(),
          "fault_ids"_a, "weight"_a, "error_probability"_a, "has_error_probability"_a, "weight_is_negative"_a, u8R"(
     Initialises a WeightedEdgeData object

     Parameters
     ----------
     fault_ids: set[int]
         A set of fault IDs
     weight: float
         The edge weight
     error_probability: float
         The probability that the edge flips. If no error_probability is associated
         with the edge, set to -1.
     has_error_probability: bool
         Whether the edge has an error_probability
     )")
     .def("__repr__", &WeightedEdgeData::repr)
     .def_readwrite("fault_ids", &WeightedEdgeData::fault_ids)
     .def_readwrite("weight", &WeightedEdgeData::weight)
     .def_readwrite("error_probability", &WeightedEdgeData::error_probability)
     .def_readwrite("has_error_probability", &WeightedEdgeData::has_error_probability);

     py::class_<MatchingGraph>(m, "MatchingGraph", u8R"(
     A matching graph to be decoded with minimum-weight perfect matching

     Examples
     --------
     >>> import math
     >>> from pymatching._cpp_mwpm import MatchingGraph, set_seed
     >>> set_seed(0)
     >>> graph = MatchingGraph()
     >>> graph.add_edge(0, 1, fault_ids={0}, weight=math.log((1-0.1)/0.1), error_probability=0.1, has_error_probability=True)
     >>> graph.add_edge(1, 2, fault_ids={1}, weight=math.log((1-0.2)/0.2), error_probability=0.2, has_error_probability=True)
     >>> graph.add_edge(2, 3, fault_ids={2}, weight=math.log((1-0.4)/0.4), error_probability=0.4, has_error_probability=True)
     >>> graph.add_edge(3, 4, fault_ids={3}, weight=math.log((1-0.05)/0.05), error_probability=0.05, has_error_probability=True)
     >>> graph.set_boundary({0, 4})
     >>> graph
     <pymatching._cpp_mwpm.MatchingGraph object with 5 nodes, 4 edges and 2 boundary nodes>
     >>> graph.get_path(1, 4)
     [1, 2, 3, 4]
     >>> graph.get_nearest_neighbours(source=2, num_neighbours=2, defect_id=[-1, 0, 1, 2, -1])
     [(3, 0.4054651081081642), (1, 1.3862943611198906)]
     )")
     .def(py::init<>(), u8R"(
     Initialises a `pymatching._cpp_mwpm.MatchingGraph`
     )")
     .def(py::init<int, std::set<int>&>(), "num_detectors"_a, "boundary"_a, u8R"(
     Initialises a `pymatching._cpp_mwpm.MatchingGraph`

     Parameters
     ----------
     num_detectors: int
         The number of detectors in the matching graph. A detector is a node that is not
         a boundary node, and has the same meaning as in Stim.
     boundary: set[int]
         The ids of the boundary nodes in the matching graph.

     )")
     .def("all_edges_have_error_probabilities",
          &MatchingGraph::AllEdgesHaveErrorProbabilities, u8R"(
          Outputs whether or not all edges have been assigned error probabilities

          Returns
          -------
          bool
              True if all edges have been assigned error probabilities, else False

          Examples
          --------
          >>> import math
          >>> from pymatching._cpp_mwpm import MatchingGraph
          >>> graph = MatchingGraph()
          >>> graph.add_edge(0, 1, {0}, math.log((1-0.1)/0.1), 0.1, True)
          >>> graph.all_edges_have_error_probabilities()
          True
          >>> graph.add_edge(1, 2, {1}, 0, -1, False)
          >>> graph.all_edges_have_error_probabilities()
          False
          )")
     .def("add_edge", &MatchingGraph::AddEdge, "node1"_a, "node2"_a, "fault_ids"_a,
          "weight"_a, "error_probability"_a=-1.0, "has_error_probability"_a=false, u8R"(
          Adds an edge to the matching graph

          Parameters
          ----------
          node1: int
              The id of the first node in the edge to be added
          node2: int
              The id of the second node in the edge to be added
          fault_ids: set[int]
              The ids of the self-inverse faults that flip if the edge is flipped
          weight: float
              The weight of the edge
          error_probability: float
              The probability that the edge is flipped. This parameter is optional
              and should be set to -1 (the default value) if no error probability
              needs to be set for the edge
          has_error_probability: bool
              Whether or not the edge has been given an error probability

          Examples
          --------
          >>> import math
          >>> from pymatching._cpp_mwpm import MatchingGraph
          >>> graph = MatchingGraph()
          >>> graph.add_edge(0, 1, {0, 1}, math.log((1-0.1)/0.1), 0.1, True)
          >>> graph.add_edge(1, 2, {2}, math.log((1-0.1)/0.1), 0.1, True)
          >>> graph.add_edge(2, 3, {3}, math.log((1-0.2)/0.2), 0.2, True)
          >>> graph.add_edge(3, 0, set(), 0, -1, False)
          >>> graph.set_boundary({0, 3})
          >>> graph.get_num_fault_ids()
          4
          >>> graph.get_num_edges()
          4
          )")
     .def("add_noise", &MatchingGraph::AddNoise, u8R"(
          Flips each edge independently with the associated error probability,
          returning the noise vector and syndrome.

          Returns
          -------
          numpy.ndarray
              A binary array (of dtype numpy.uint8) specifying whether each
              fault has occurred. Element i is one if the fault with
              fault ID `i`` has been flipped, and is zero otherwise.
          numpy.ndarray
              A binary array (of dtype numpy.uint8) with length equal to the
              number of nodes in the matching graph, specifying the syndrome.
              Element i is one if node i is a defect, and zero otherwise. Note
              that boundary nodes are never defects (their syndrome is zero).

          Examples
          --------
          >>> import math
          >>> from pymatching._cpp_mwpm import MatchingGraph, set_seed
          >>> set_seed(0)
          >>> graph = MatchingGraph()
          >>> graph.add_edge(0, 1, {0}, math.log((1-0.4)/0.4), 0.4, True)
          >>> graph.add_edge(1, 2, {1}, math.log((1-0.0001)/0.0001), 0.0001, True)
          >>> graph.add_edge(2, 3, {2}, math.log((1-0.45)/0.45), 0.45, True)
          >>> graph.add_noise()
          (array([0, 0, 0], dtype=uint8), array([0, 0, 0, 0], dtype=uint8))
          >>> graph.add_noise()
          (array([0, 0, 1], dtype=uint8), array([0, 0, 1, 1], dtype=uint8))
          >>> graph.add_noise()
          (array([1, 0, 1], dtype=uint8), array([1, 1, 1, 1], dtype=uint8))
          >>> graph.add_noise()
          (array([0, 0, 0], dtype=uint8), array([0, 0, 0, 0], dtype=uint8))
     )")
     .def("get_boundary", &MatchingGraph::GetBoundary, u8R"(
     Get the ids of the boundary nodes

     Returns
     -------
     set[int]
         The ids of the boundary nodes

     Examples
     --------
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.set_boundary({1,3,5})
     >>> graph.get_boundary()
     {1, 3, 5}
     )")
     .def("set_boundary", &MatchingGraph::SetBoundary, "boundary"_a, u8R"(
     Set the ids of the boundary nodes

     Parameters
     -------
     boundary: set[int]
         The ids of the boundary nodes

     Examples
     --------
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.set_boundary({2,3,4})
     >>> graph.get_boundary()
     {2, 3, 4}
     )")
     .def("get_edges", &MatchingGraph::GetEdges, u8R"(
     Get the edges and edge data in the MatchingGraph

     Returns
     -------
     list[tuple[int, int, WeightedEdgeData]]
         A list of edges. Each edges is a tuple (n1, n2, edge_data) where n1 and n2
         are the ids of the first and second nodes in the edge, respectively, and
         edge_data is the WeightedEdgeData associated with the edge.

     Examples
     --------
     >>> import math
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.add_edge(0, 1, {0}, math.log((1-0.4)/0.4), 0.4, True)
     >>> graph.add_edge(1, 2, {1}, math.log((1-0.3)/0.3), 0.3, True)
     >>> graph.get_edges()
     [(0, 1, pymatching._cpp_mwpm.WeightedEdgeData({0}, 0.405465, 0.4, 1)), (1, 2, pymatching._cpp_mwpm.WeightedEdgeData({1}, 0.847298, 0.3, 1))]
     )")
     .def("get_nearest_neighbours", &MatchingGraph::GetNearestNeighbours,
          "source"_a, "num_neighbours"_a, "defect_id"_a, u8R"(
          Find the nearest `num_neighbours` defects from a source node in the matching graph,
          using a modified Dijkstra algorithm (local Dijkstra).

          Parameters
          ----------
          source: int
              The index of the source node
          num_neighbours: int
              The maximum number of defects to find in the matching graph (excluding the source node)
          defect_id: list[int]
              A list of length `MatchingGraph.get_num_nodes()` specifying the defect id of each node
              in the matching graph. Node `i` satisfies `defect_id[i]>=0` if it is a defect and
              `defect_id[i]==-1` otherwise. The value of `defect_id[i]` (if not -1) is the index
              of the corresponding defect in the syndrome vector.

          Returns
          -------
          list[tuple[int, float]]
              A list of tuples of the form `(i, d)` where `i` is a node id of
              a defect, and `d` is its distance from the source node (the sum of the
              weights of the edges along the shortest path from the source to node `i`).
              The list has `num_neighbours` elements/tuples, and the nodes in the list
              are the `num_neighbours` closest defects from the source node.

          Examples
          --------
          >>> from pymatching._cpp_mwpm import MatchingGraph
          >>> graph = MatchingGraph()
          >>> graph.add_edge(0, 1, {0}, 0.1)
          >>> graph.add_edge(1, 2, {1}, 0.2)
          >>> graph.add_edge(2, 3, {2}, 0.1)
          >>> graph.add_edge(3, 4, {3}, 0.1)
          >>> graph.add_edge(4, 5, {4}, 0.1)
          >>> graph.get_nearest_neighbours(2, 0, [-1, 0, 1, 2, -1, 3])
          []
          >>> graph.get_nearest_neighbours(2, 1, [-1, 0, 1, 2, -1, 3])
          [(3, 0.1)]
          >>> graph.get_nearest_neighbours(2, 2, [-1, 0, 1, 2, -1, 3])
          [(3, 0.1), (1, 0.2)]
          >>> graph.get_nearest_neighbours(2, 3, [-1, 0, 1, 2, -1, 3])
          [(3, 0.1), (1, 0.2), (5, 0.30000000000000004)]
          )")
     .def("get_path", &MatchingGraph::GetPath, "source"_a, "target"_a, u8R"(
          Find the nodes along the shortest path from a source to a target
          node using Dijkstra's algorithm.

          Parameters
          ----------
          source: int
              The id of the source vertex
          target: int
              The id of the target vertex

          Returns
          -------
          list[int]
              A list of ids of the nodes along the shortest path from source to
              target (including the source and target nodes). Elements `i` and `i+1`
              of the list are nodes in the `i`th edge along the shortest path from
              source to target.

          Examples
          --------
          >>> from pymatching._cpp_mwpm import MatchingGraph
          >>> graph = MatchingGraph()
          >>> graph.add_edge(0, 1, {0}, 1)
          >>> graph.add_edge(1, 2, {1}, 1)
          >>> graph.add_edge(2, 3, {2}, 1)
          >>> graph.add_edge(3, 4, {3}, 1)
          >>> graph.add_edge(4, 5, {4}, 1)
          >>> graph.get_path(2, 5)
          [2, 3, 4, 5]
          >>> graph.get_path(4, 1)
          [4, 3, 2, 1]
          >>> graph.get_path(5, 0)
          [5, 4, 3, 2, 1, 0]
     )")
     .def("distance", &MatchingGraph::Distance, "node1"_a, "node2"_a, u8R"(
     Get the distance between node1 and node2 using the precomputed all-pairs shortest paths
     computed using Dijkstra. If the all-pairs-shortest-paths have not yet been computed, this
     function will also compute these.

     Parameters
     ----------
     node1: int
         The id of the first node
     node2: int
         The id of the second node

     Returns
     -------
     float
         The sum of the weights of the edges along the shortest path from node1 to node2

     Examples
     --------
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.add_edge(0, 1, {0}, 0.12)
     >>> graph.add_edge(1, 2, {1}, 0.11)
     >>> graph.add_edge(2, 3, {2}, 0.3)
     >>> graph.distance(0, 3)
     0.53
     >>> graph.distance(1, 2)
     0.11
     >>> graph.distance(3, 2)
     0.3
     )")
     .def("shortest_path", &MatchingGraph::ShortestPath, "node1"_a, "node2"_a, u8R"(
     Find the shortest path between node1 and node2, using the precomputed all-pairs shortest paths
     computed using Dijkstra's algorithm. If the all-pairs-shortest-paths have not yet been computed, this
     function will also compute these.

     Parameters
     ----------
     node1: int
         The id of the first node
     node2: int
         The id of the second node

     Returns
     -------
     list[int]
         A list of ids of the nodes along the shortest path from source to
         target (including the source and target nodes). Elements `i` and `i+1`
         of the list are nodes in the `i`th edge along the shortest path from
         source to target.

     Examples
     --------
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.add_edge(0, 1, {0}, 1)
     >>> graph.add_edge(1, 2, {1}, 1)
     >>> graph.add_edge(2, 3, {2}, 1)
     >>> graph.add_edge(3, 4, {3}, 1)
     >>> graph.shortest_path(0, 4)
     [0, 1, 2, 3, 4]
     >>> graph.shortest_path(3, 2)
     [3, 2]
     >>> graph.shortest_path(1, 2)
     [1, 2]
     )")
     .def("fault_ids", &MatchingGraph::FaultIDs, "node1"_a, "node2"_a, u8R"(
     Returns the fault_ids associated with the edge (node1, node2)

     Parameters
     ----------
     node1: int
         The id of the first node
     node2: int
         The id of the second node

     Returns
     -------
     set[int]
         The fault_ids associated with the edge (node1, node2)

     Examples
     --------
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.add_edge(0, 1, {0}, 1)
     >>> graph.add_edge(1, 2, {1, 2, 3}, 1)
     >>> graph.fault_ids(0, 1)
     {0}
     >>> graph.fault_ids(1, 2)
     {1, 2, 3}
     >>> graph.fault_ids(1, 0)
     {0}
     )")
     .def("get_num_fault_ids", &MatchingGraph::GetNumFaultIDs, u8R"(
     Returns the number of distinct fault_ids associated with edges in the matching graph.

     Returns
     -------
     int
         The number of distinct fault_ids in the matching graph

     Examples
     --------
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.add_edge(0, 1, {0}, 1)
     >>> graph.add_edge(1, 2, {1, 2, 3}, 1)
     >>> graph.add_edge(2, 3, {4, 5}, 1)
     >>> graph.get_num_fault_ids()
     6
     )")
     .def("get_num_nodes", &MatchingGraph::GetNumNodes, u8R"(
     Returns the number of nodes in the matching graph

     Returns
     -------
     int
         The number of nodes in the matching graph

     Examples
     --------
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.add_edge(0, 1, {0}, 1)
     >>> graph.add_edge(1, 2, {1}, 1)
     >>> graph.add_edge(2, 3, {2}, 1)
     >>> graph.get_num_nodes()
     4
     )")
     .def("get_num_edges", &MatchingGraph::GetNumEdges, u8R"(
     Get the number of edges in the matching graph

     Returns
     -------
     int
         The number of edges in the matching graph

     Examples
     --------
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.add_edge(0, 1, {0}, 1)
     >>> graph.add_edge(1, 2, {1}, 1)
     >>> graph.add_edge(2, 3, {2}, 1)
     >>> graph.get_num_edges()
     3
     )")
     .def("compute_all_pairs_shortest_paths", &MatchingGraph::ComputeAllPairsShortestPaths, u8R"(
     Computes the shortest paths between all pairs of nodes in the matching graph using Dijkstra's
     algorithm. Note that this method is very memory intensive and is not used for local matching,
     only for exact matching.
     )")
     .def("has_computed_all_pairs_shortest_paths", &MatchingGraph::HasComputedAllPairsShortestPaths, u8R"(
     Returns whether or not the all-pairs shortest paths have already been computed.

     Examples
     --------
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.add_edge(0, 1, {0}, 0.5)
     >>> graph.add_edge(1, 2, {1}, 1)
     >>> graph.has_computed_all_pairs_shortest_paths()
     False
     >>> graph.compute_all_pairs_shortest_paths()
     >>> graph.has_computed_all_pairs_shortest_paths()
     True
     >>> graph.add_edge(2, 3, {2}, 0.8)
     >>> graph.has_computed_all_pairs_shortest_paths()
     False
     >>> graph.compute_all_pairs_shortest_paths()
     >>> graph.has_computed_all_pairs_shortest_paths()
     True
     )")
     .def("get_num_connected_components", &MatchingGraph::GetNumConnectedComponents, u8R"(
     Get the number of connected components in the matching graph

     Returns
     -------
     int
         The number of connected components in the matching graph

     Examples
     --------
     >>> from pymatching._cpp_mwpm import MatchingGraph
     >>> graph = MatchingGraph()
     >>> graph.add_edge(0, 1, {0}, 1)
     >>> graph.add_edge(1, 2, {1}, 1)
     >>> graph.get_num_connected_components()
     1
     >>> graph.add_edge(3, 4, {2}, 1)
     >>> graph.get_num_connected_components()
     2
     >>> graph.add_edge(5, 6, {3}, 1)
     >>> graph.add_edge(6, 7, {4}, 1)
     >>> graph.get_num_connected_components()
     3
     )")
     .def("__repr__", &MatchingGraph::repr);

     m.def("randomize", &randomize, u8R"(
        Choose a random seed using std::random_device

        Examples
        --------
            >>> import pymatching
            >>> pymatching.randomize()
     )");
     m.def("set_seed", &set_seed, "seed"_a, u8R"(
        Sets the seed of the random number generator

        Parameters
        ----------
        seed: int
            The seed for the random number generator (must be non-negative)

        Examples
        --------
        >>> import pymatching
        >>> pymatching.set_seed(10)

     )");
     m.def("rand_float", &rand_float, "from"_a, "to"_a, u8R"(
        Generate a floating point number chosen uniformly at random
        over the interval between `from` and `to`

        Parameters
        ----------
        from: float
            Smallest float that can be drawn from the distribution

        to: float
            Largest float that can be drawn from the distribution

        Returns
        -------
        float
            The random float
     )");

     py::class_<MatchingResult>(m, "MatchingResult")
     .def(py::init<>())
     .def(py::init<py::array_t<std::uint8_t>, double>(), "correction"_a, "weight"_a)
     .def_readwrite("correction", &MatchingResult::correction)
     .def_readwrite("weight", &MatchingResult::weight)
     .def("__repr__", &MatchingResult::repr);

     py::register_exception<BlossomFailureException>(m, "BlossomFailureException", PyExc_RuntimeError);

     m.def("local_matching", &LocalMatching,
           "graph"_a, "defects"_a, "num_neighbours"_a=30,
           "return_weight"_a=false, "max_attempts"_a=10, u8R"(
            Decode using local matching.

            Parameters
            ----------
            graph: MatchingGraph
                The matching graph to be used to decode the syndrome
            defects: np.ndarray[int]
                A numpy array of integers giving the indices of the -1 measurements
                (defects) in the syndrome. i.e. This is an array of IDs of nodes with
                a non-trivial syndrome (detectors that have fired).
            num_neighbours: int
                Number of closest neighbours (with non-trivial syndrome) of each matching
                graph node to consider when decoding. `num_neighbours` corresponds to
                the parameter `m` in the local matching algorithm in the paper:
                https://arxiv.org/abs/2105.13082 is used, and `num_neighbours`
                It is recommended to set `num_neighbours` to at least 20 for
                decoding performance to closely match that of exact matching.
            return_weight: bool
                If True, also return the weight of the matching found. By default False
            max_attempts: int
                The blossom algorithm can very occasionally fail to find a solution if a
                perfect matching does not exist in the graph derived from nodes only with
                non-trivial syndromes (called a syndrome graph in https://arxiv.org/abs/2105.13082),
                since the syndrome graph is not a complete graph in local matching.
                If this happens, `num_neighbours` is doubled `max_attempts` times until
                a solution is found. It is highly unlikely that more than one attempt is
                required, but if a solution is not found after `max_attempts` tries
                (doubling `num_neighbours` each time), then a `pymatching.BlossomFailureException`
                is raised. By default 10

            Returns
            -------
            MatchingResult
                The recovery operator (and, optionally, weight of the matching).
                `MatchingResult.correction[i]` is 1 if fault_ids 1
                should be flipped when applying the minimum weight correction, and 0
                otherwise. `MatchingResult.weight` gives the sum of the weights of the
                edges included in the minimum-weight perfect matching correction if
                `return_weight=True`, and is -1 otherwise.

            Examples
            --------
            >>> from pymatching._cpp_mwpm import local_matching, MatchingGraph
            >>> graph = MatchingGraph()
            >>> graph.add_edge(0, 1, {0}, 1.0)
            >>> graph.add_edge(1, 2, {1}, 1.0)
            >>> graph.add_edge(2, 0, {2}, 1.0)
            >>> res = local_matching(graph, [0, 1])
            >>> res.correction
            array([1, 0, 0], dtype=uint8)

            By setting `return_weight=True`, the weight of the matching is also
            returned:
            >>> from pymatching._cpp_mwpm import local_matching, MatchingGraph
            >>> graph = MatchingGraph()
            >>> graph.add_edge(0, 1, {0}, 0.23)
            >>> graph.add_edge(1, 2, {1}, 0.34)
            >>> graph.add_edge(2, 3, {2}, 0.91)
            >>> graph.add_edge(3, 0, {3}, 0.86)
            >>> res = local_matching(graph, [0, 2], return_weight=True)
            >>> res.correction
            array([1, 1, 0, 0], dtype=uint8)
            >>> res.weight
            0.5700000000000001
           )");
     m.def("exact_matching", &ExactMatching, "graph"_a, "defects"_a,
           "return_weight"_a=false, u8R"(
            Decode using exact matching

            Parameters
            ----------
            graph: MatchingGraph
                The matching graph to be used to decode the syndrome
            defects: np.ndarray[int]
                A numpy array of integers giving the indices of the -1 measurements
                (defects) in the syndrome. i.e. This is an array of IDs of nodes with
                a non-trivial syndrome (detectors that have fired).
            return_weight: bool
                If True, also return the weight of the matching found. By default False

            Returns
            -------
            MatchingResult
                The recovery operator (and, optionally, weight of the matching).
                `MatchingResult.correction[i]` is 1 if fault_ids 1
                should be flipped when applying the minimum weight correction, and 0
                otherwise. `MatchingResult.weight` gives the sum of the weights of the
                edges included in the minimum-weight perfect matching correction if
                `return_weight=True`, and is -1 otherwise.

            Examples
            --------
            >>> from pymatching._cpp_mwpm import exact_matching, MatchingGraph
            >>> graph = MatchingGraph()
            >>> graph.add_edge(0, 1, {0}, 0.8)
            >>> graph.add_edge(1, 2, {1}, 0.5)
            >>> graph.add_edge(2, 3, {2, 3}, 0.9)
            >>> graph.add_edge(3, 0, {4}, 1.1)
            >>> res = exact_matching(graph, [1, 3])
            >>> res.correction
            array([0, 1, 1, 1, 0], dtype=uint8)
           )");
}
