#ifndef PYMATCHING2_USER_GRAPH_PYBIND_H
#define PYMATCHING2_USER_GRAPH_PYBIND_H

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pymatching/fill_match/driver/user_graph.h"

namespace py = pybind11;

namespace pm_pybind {

py::class_<pm::UserGraph> pybind_user_graph(py::module &m);

void pybind_user_graph_methods(py::module &m, py::class_<pm::UserGraph> &g);

}  // namespace pm_pybind

#endif  // PYMATCHING2_USER_GRAPH_PYBIND_H
