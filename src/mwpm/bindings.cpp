#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include"mwpm.h"
#include <PerfectMatching.h>

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

    py::class_<GraphData>(m, "GraphData")
    .def(py::init<>())
    .def_readwrite("g", &GraphData::g)
    .def_readwrite("qubit", &GraphData::qubit);

    m.def("_stabiliser_graph", &StabiliserGraph, "indices"_a, "num_stabilisers"_a);

    py::class_<BFSResult>(m, "BFSResult")
    .def(py::init<>())
    .def_readwrite("distance", &BFSResult::distance)
    .def_readwrite("parent", &BFSResult::parent);

    py::class_<APSPResult>(m, "APSPResult")
    .def(py::init<>())
    .def_readwrite("distances", &APSPResult::distances)
    .def_readwrite("parents", &APSPResult::parents);

    py::class_<StabiliserGraphCls>(m, "StabiliserGraph")
    .def(py::init<const py::array_t<int>&, int, int>(), "indices"_a, "num_stabilisers"_a, "num_qubits"_a)
    .def_readwrite("adj_list", &StabiliserGraphCls::adj_list)
    .def_readwrite("qubit_ids", &StabiliserGraphCls::qubit_ids)
    .def_readwrite("shortest_paths", &StabiliserGraphCls::shortest_paths)
    .def_readwrite("num_qubits", &StabiliserGraphCls::num_qubits)
    .def_readwrite("num_stabilisers", &StabiliserGraphCls::num_stabilisers)
    .def("distance", &StabiliserGraphCls::Distance)
    .def("shortest_path", &StabiliserGraphCls::ShortestPath)
    .def("qubit_id", &StabiliserGraphCls::QubitID);

    m.def("breadth_first_search", &BreadthFirstSearch, "g"_a, "source"_a);
    m.def("all_pairs_shortest_path", &AllPairsShortestPath, "g"_a);
    m.def("shortest_path", &GetShortestPath, "parent"_a, "dest"_a);
    m.def("_decode", &Decode, "sg"_a, "defects"_a);
}
