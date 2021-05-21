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
#include "stabiliser_graph.h"
#include "weighted_stabiliser_graph.h"
#include "rand_gen.h"
#include "lemon_mwpm.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(_cpp_mwpm, m) {
     py::class_<IStabiliserGraph>(m, "IStabiliserGraph")
     .def("distance", &IStabiliserGraph::Distance, "node1"_a, "node2"_a)
     .def("space_time_distance", &IStabiliserGraph::SpaceTimeDistance, "node1"_a, "node2"_a)
     .def("shortest_path", &IStabiliserGraph::ShortestPath, "node1"_a, "node2"_a)
     .def("space_time_shortest_path", &IStabiliserGraph::SpaceTimeShortestPath, "node1"_a, "node2"_a)
     .def("qubit_ids", &IStabiliserGraph::QubitIDs, "node1"_a, "node2"_a)
     .def("get_num_qubits", &IStabiliserGraph::GetNumQubits)
     .def("get_num_nodes", &IStabiliserGraph::GetNumNodes)
     .def("get_num_edges", &IStabiliserGraph::GetNumEdges)
     .def("compute_all_pairs_shortest_paths", &IStabiliserGraph::ComputeAllPairsShortestPaths)
     .def("has_computed_all_pairs_shortest_paths", &IStabiliserGraph::HasComputedAllPairsShortestPaths)
     .def("get_num_connected_components", &IStabiliserGraph::GetNumConnectedComponents);

     py::class_<WeightedEdgeData>(m, "WeightedEdgeData")
     .def_readwrite("qubit_ids", &WeightedEdgeData::qubit_ids)
     .def_readwrite("weight", &WeightedEdgeData::weight)
     .def_readwrite("error_probability", &WeightedEdgeData::error_probability)
     .def_readwrite("has_error_probability", &WeightedEdgeData::has_error_probability);

     py::class_<WeightedStabiliserGraph, IStabiliserGraph>(m, "WeightedStabiliserGraph")
     .def(py::init<int, std::vector<int>&>(), "num_stabilisers"_a, "boundary"_a)
     .def_readwrite("all_predecessors", &WeightedStabiliserGraph::all_predecessors)
     .def_readwrite("all_distances", &WeightedStabiliserGraph::all_distances)
     .def_readwrite("all_edges_have_error_probabilities", 
                    &WeightedStabiliserGraph::all_edges_have_error_probabilities)
     .def_readwrite("_distances", &WeightedStabiliserGraph::_distances)
     .def_readwrite("_predecessors", &WeightedStabiliserGraph::_predecessors)
     .def("add_edge", &WeightedStabiliserGraph::AddEdge, "node1"_a, "node2"_a, "qubit_ids"_a, 
          "weight"_a, "error_probability"_a=-1.0, "has_error_probability"_a=false)
     .def("add_noise", &WeightedStabiliserGraph::AddNoise)
     .def("get_boundary", &WeightedStabiliserGraph::GetBoundary)
     .def("set_boundary", &WeightedStabiliserGraph::SetBoundary, "boundary"_a)
     .def("get_edges", &WeightedStabiliserGraph::GetEdges)
     .def("get_nearest_neighbours", &WeightedStabiliserGraph::GetNearestNeighbours,
          "source"_a, "num_neighbours"_a, "defect_id"_a)
     .def("get_path", &WeightedStabiliserGraph::GetPath, "source"_a, "target"_a);

     m.def("randomize", &randomize);
     m.def("set_seed", &set_seed, "s"_a);
     m.def("rand_float", &rand_float, "from"_a, "to"_a);

     py::class_<MatchingResult>(m, "MatchingResult")
     .def_readwrite("correction", &MatchingResult::correction)
     .def_readwrite("weight", &MatchingResult::weight);

     m.def("decode_match_neighbourhood", &LemonDecodeMatchNeighbourhood,
           "sg"_a, "defects"_a, "num_neighbours"_a=20,
           "return_weight"_a=false);
     m.def("decode", &LemonDecode, "sg"_a, "defects"_a,
           "return_weight"_a=false);
}
