#include "pymatching/fill_match/driver/user_graph.pybind.h"
#include "pybind11/pybind11.h"

namespace py = pybind11;
using namespace py::literals;



PYBIND11_MODULE(PYMATCHING_PYBIND11_MODULE_NAME, m) {
    auto matching_graph = pm_pybind::pybind_user_graph(m);
    pm_pybind::pybind_user_graph_methods(m, matching_graph);
}