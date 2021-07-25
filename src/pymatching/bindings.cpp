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
     .def_readwrite("qubit_ids", &WeightedEdgeData::qubit_ids)
     .def_readwrite("weight", &WeightedEdgeData::weight)
     .def_readwrite("error_probability", &WeightedEdgeData::error_probability)
     .def_readwrite("has_error_probability", &WeightedEdgeData::has_error_probability);

     py::class_<MatchingGraph>(m, "MatchingGraph")
     .def(py::init<>())
     .def(py::init<int, std::set<int>&>(), "num_detectors"_a, "boundary"_a)
     .def("all_edges_have_error_probabilities",
          &MatchingGraph::AllEdgesHaveErrorProbabilities)
     .def("add_edge", &MatchingGraph::AddEdge, "node1"_a, "node2"_a, "qubit_ids"_a,
          "weight"_a, "error_probability"_a=-1.0, "has_error_probability"_a=false)
     .def("add_noise", &MatchingGraph::AddNoise)
     .def("get_boundary", &MatchingGraph::GetBoundary)
     .def("set_boundary", &MatchingGraph::SetBoundary, "boundary"_a)
     .def("get_edges", &MatchingGraph::GetEdges)
     .def("get_nearest_neighbours", &MatchingGraph::GetNearestNeighbours,
          "source"_a, "num_neighbours"_a, "defect_id"_a)
     .def("get_path", &MatchingGraph::GetPath, "source"_a, "target"_a)
     .def("distance", &MatchingGraph::Distance, "node1"_a, "node2"_a)
     .def("shortest_path", &MatchingGraph::ShortestPath, "node1"_a, "node2"_a)
     .def("qubit_ids", &MatchingGraph::QubitIDs, "node1"_a, "node2"_a)
     .def("get_num_qubits", &MatchingGraph::GetNumQubits)
     .def("get_num_nodes", &MatchingGraph::GetNumNodes)
     .def("get_num_edges", &MatchingGraph::GetNumEdges)
     .def("compute_all_pairs_shortest_paths", &MatchingGraph::ComputeAllPairsShortestPaths)
     .def("has_computed_all_pairs_shortest_paths", &MatchingGraph::HasComputedAllPairsShortestPaths)
     .def("get_num_connected_components", &MatchingGraph::GetNumConnectedComponents);

     m.def("randomize", &randomize);
     m.def("set_seed", &set_seed, "s"_a);
     m.def("rand_float", &rand_float, "from"_a, "to"_a);

     py::class_<MatchingResult>(m, "MatchingResult")
     .def_readwrite("correction", &MatchingResult::correction)
     .def_readwrite("weight", &MatchingResult::weight);

     py::register_exception<BlossomFailureException>(m, "BlossomFailureException", PyExc_RuntimeError);

     m.def("local_matching", &LocalMatching,
           "sg"_a, "defects"_a, "num_neighbours"_a=30,
           "return_weight"_a=false, "max_attempts"_a=10);
     m.def("exact_matching", &LemonDecode, "sg"_a, "defects"_a,
           "return_weight"_a=false);
}
