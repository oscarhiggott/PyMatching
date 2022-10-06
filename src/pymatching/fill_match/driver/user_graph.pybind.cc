
#include "pymatching/fill_match/driver/user_graph.pybind.h"

#include <pybind11/numpy.h>

#include "pybind11/pybind11.h"
#include "pymatching/fill_match/driver/mwpm_decoding.h"

using namespace py::literals;

py::class_<pm::UserGraph> pm_pybind::pybind_user_graph(py::module &m) {
    auto g = py::class_<pm::UserGraph>(m, "MatchingGraph");
    return g;
}

void pm_pybind::pybind_user_graph_methods(py::module &m, py::class_<pm::UserGraph> &g) {
    g.def(py::init<>());
    g.def(py::init<size_t>(), "num_nodes"_a);
    g.def(
        "add_edge",
        &pm::UserGraph::add_or_merge_edge, "node1"_a, "node2"_a, "observables"_a, "weight"_a, "error_probability"_a);
    g.def(
        "add_boundary_edge",
        &pm::UserGraph::add_or_merge_boundary_edge,
        "node"_a,
        "observables"_a,
        "weight"_a,
        "error_probability"_a);
    g.def("set_boundary", &pm::UserGraph::set_boundary, "boundary"_a);
    g.def("get_boundary", &pm::UserGraph::get_boundary);
    g.def("get_num_observables", &pm::UserGraph::get_num_observables);
    g.def("get_num_nodes", &pm::UserGraph::get_num_nodes);
    g.def("get_num_edges", &pm::UserGraph::get_num_edges);
    g.def("get_num_detectors", &pm::UserGraph::get_num_detectors);
    g.def("all_edges_have_error_probabilities", &pm::UserGraph::all_edges_have_error_probabilities);
    g.def("add_noise", [](pm::UserGraph &self) {
        auto error_vec = new std::vector<uint8_t>(self.get_num_observables(), 0);
        auto syndrome_vec = new std::vector<uint8_t>(self.get_num_nodes(), 0);
        self.add_noise(error_vec->data(), syndrome_vec->data());

        auto syndrome_capsule = py::capsule(syndrome_vec, [](void *syndrome) {
            delete reinterpret_cast<std::vector<uint8_t> *>(syndrome);
        });
        py::array_t<int> syndrome_arr =
            py::array_t<uint8_t>(syndrome_vec->size(), syndrome_vec->data(), syndrome_capsule);
        auto err_capsule = py::capsule(error_vec, [](void *error) {
            delete reinterpret_cast<std::vector<uint8_t> *>(error);
        });
        py::array_t<int> error_arr = py::array_t<uint8_t>(error_vec->size(), error_vec->data(), err_capsule);

        std::pair<py::array_t<std::uint8_t>, py::array_t<std::uint8_t>> res = {error_arr, syndrome_arr};
        return res;
    });
    g.def("decode", [](pm::UserGraph &self, const py::array_t<uint64_t> &detection_events) {
        std::vector<uint64_t> detection_events_vec(
            detection_events.data(), detection_events.data() + detection_events.size());
        auto &mwpm = self.get_mwpm();
        auto obs_crossed = new std::vector<uint8_t>(self.get_num_observables(), 0);
        pm::total_weight_int weight = 0;
        pm::decode_detection_events(mwpm, detection_events_vec, obs_crossed->data(), weight);
        double rescaled_weight = (double)weight / mwpm.flooder.graph.normalising_constant;

        auto err_capsule = py::capsule(obs_crossed, [](void *x) {
            delete reinterpret_cast<std::vector<uint8_t> *>(x);
        });
        py::array_t<uint8_t> obs_crossed_arr =
            py::array_t<uint8_t>(obs_crossed->size(), obs_crossed->data(), err_capsule);
        std::pair<py::array_t<std::uint8_t>, double> res = {obs_crossed_arr, rescaled_weight};
        return res;
    });
    g.def("get_edges", &pm::UserGraph::get_edges);
}