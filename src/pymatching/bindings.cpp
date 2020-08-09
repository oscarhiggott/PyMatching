#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include"mwpm.h"
#include <PerfectMatching.h>
#include "stabiliser_graph.h"
#include "weighted_stabiliser_graph.h"
#include "rand_gen.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(_cpp_mwpm, m) {
    py::class_<PerfectMatching>(m, "PerfectMatching")
    .def(py::init<int, int>(), "nodeNum"_a, "edgeNumMax"_a)
    .def("add_edge", &PerfectMatching::AddEdge, 
         "i"_a, "j"_a, "cost"_a)
    .def("solve", &PerfectMatching::Solve, "finish"_a=true)
    .def("get_match", &PerfectMatching::GetMatch, "i"_a)
    .def("get_solution", &PerfectMatching::GetSolution, "e"_a)
    .def_readwrite("options", &PerfectMatching::options);

    py::class_<PerfectMatching::Options>(m, "Options")
    .def(py::init<>())
    .def_readwrite("fractional_jumpstart", &PerfectMatching::Options::fractional_jumpstart)
    .def_readwrite("dual_greedy_update_option", &PerfectMatching::Options::dual_greedy_update_option)
    .def_readwrite("dual_LP_threshold", &PerfectMatching::Options::dual_LP_threshold)
    .def_readwrite("update_duals_before", &PerfectMatching::Options::update_duals_before)
    .def_readwrite("update_duals_after", &PerfectMatching::Options::update_duals_after)
    .def_readwrite("single_tree_threshold", &PerfectMatching::Options::single_tree_threshold)
    .def_readwrite("verbose", &PerfectMatching::Options::verbose);

    py::class_<IStabiliserGraph>(m, "IStabiliserGraph")
    .def("distance", &IStabiliserGraph::Distance, "node1"_a, "node2"_a)
    .def("space_time_distance", &IStabiliserGraph::SpaceTimeDistance, "node1"_a, "node2"_a)
    .def("shortest_path", &IStabiliserGraph::ShortestPath, "node1"_a, "node2"_a)
    .def("space_time_shortest_path", &IStabiliserGraph::SpaceTimeShortestPath, "node1"_a, "node2"_a)
    .def("qubit_id", &IStabiliserGraph::QubitID, "node1"_a, "node2"_a)
    .def("get_num_qubits", &IStabiliserGraph::GetNumQubits)
    .def("get_num_stabilisers", &IStabiliserGraph::GetNumStabilisers)
    .def("compute_all_pairs_shortest_paths", &IStabiliserGraph::ComputeAllPairsShortestPaths);

    py::class_<WeightedStabiliserGraph, IStabiliserGraph>(m, "WeightedStabiliserGraph")
    .def(py::init<int, std::vector<int>&>(), "num_stabilisers"_a, "boundary"_a)
    .def(py::init<const py::array_t<int>&, const py::array_t<double>&, std::vector<int>&>(), 
        "indices"_a, "weights"_a, "boundary"_a)
    .def(py::init<const py::array_t<int>&, const py::array_t<double>&, 
            const py::array_t<double>&, std::vector<int>&>(), 
            "indices"_a, "weights"_a, "error_probabilities"_a, "boundary"_a)
    .def_readwrite("all_predecessors", &WeightedStabiliserGraph::all_predecessors)
    .def_readwrite("all_distances", &WeightedStabiliserGraph::all_distances)
    .def_readwrite("all_edges_have_error_probabilities", 
                    &WeightedStabiliserGraph::all_edges_have_error_probabilities)
    .def("add_edge", &WeightedStabiliserGraph::AddEdge, "node1"_a, "node2"_a, "qubit_id"_a, 
         "weight"_a, "error_probability"_a=-1.0, "has_error_probability"_a=false)
    .def("add_noise", &WeightedStabiliserGraph::AddNoise)
    .def("get_boundary", &WeightedStabiliserGraph::GetBoundary)
    .def("set_boundary", &WeightedStabiliserGraph::SetBoundary, "boundary"_a)
    .def("get_nearest_neighbours", &WeightedStabiliserGraph::GetNearestNeighbours,
         "source"_a, "num_neighbours"_a, "defect_id"_a)
    .def("get_path", &WeightedStabiliserGraph::GetPath, "source"_a, "target"_a);

    m.def("decode", &Decode, "sg"_a, "defects"_a);
    m.def("decode_match_neighbourhood", &DecodeMatchNeighbourhood,
          "sg"_a ,"defects"_a, "num_neighbours"_a);
    m.def("randomize", &randomize);
    m.def("set_seed", &set_seed, "s"_a);
    m.def("rand_float", &rand_float, "from"_a, "to"_a);
}
