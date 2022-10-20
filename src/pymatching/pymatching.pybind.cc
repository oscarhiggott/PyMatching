#include "pybind11/pybind11.h"
#include "pymatching/sparse_blossom/driver/namespaced_main.h"
#include "pymatching/sparse_blossom/driver/user_graph.pybind.h"
#include "pymatching/rand/rand_gen.pybind.h"

namespace py = pybind11;
using namespace py::literals;

int pymatching_main(const std::vector<std::string> &args) {
    std::vector<const char *> argv;
    argv.push_back("pymatching.main");
    for (const auto &arg : args) {
        argv.push_back(arg.data());
    }
    return pm::main(argv.size(), argv.data());
}

PYBIND11_MODULE(_cpp_pymatching, m) {
    auto matching_graph = pm_pybind::pybind_user_graph(m);
    pm_pybind::pybind_user_graph_methods(m, matching_graph);
    pm_pybind::pybind_rand_gen_methods(m);
    m.def("main", &pymatching_main, pybind11::kw_only(), pybind11::arg("command_line_args"), R"pbdoc(
Runs the command line tool version of pymatching with the given arguments.
)pbdoc");
}