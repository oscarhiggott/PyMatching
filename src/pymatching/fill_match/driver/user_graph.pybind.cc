
#include "pymatching/fill_match/driver/user_graph.pybind.h"

#include <pybind11/numpy.h>

#include "pybind11/pybind11.h"
#include "pymatching/fill_match/driver/mwpm_decoding.h"
#include "stim.h"

using namespace py::literals;

py::class_<pm::UserGraph> pm_pybind::pybind_user_graph(py::module &m) {
    auto g = py::class_<pm::UserGraph>(m, "MatchingGraph");
    return g;
}

pm::MERGE_STRATEGY merge_strategy_from_string(const std::string &merge_strategy) {
    static std::unordered_map<std::string, pm::MERGE_STRATEGY> const table = {
        {"disallow", pm::DISALLOW},
        {"independent", pm::INDEPENDENT},
        {"smallest-weight", pm::SMALLEST_WEIGHT},
        {"keep-original", pm::KEEP_ORIGINAL},
        {"replace", pm::REPLACE}};
    auto it = table.find(merge_strategy);
    if (it != table.end()) {
        return it->second;
    } else {
        throw std::invalid_argument("Merge strategy \"" + merge_strategy + "\" not recognised.");
    }
}

void pm_pybind::pybind_user_graph_methods(py::module &m, py::class_<pm::UserGraph> &g) {
    g.def(py::init<>());
    g.def(py::init<size_t>(), "num_nodes"_a);
    g.def(py::init<size_t, size_t>(), "num_nodes"_a, "num_fault_ids"_a);
    g.def(
        "add_edge",
        [](pm::UserGraph &self,
           int64_t node1,
           int64_t node2,
           const std::set<size_t> &observables,
           double weight,
           double error_probability,
           const std::string &merge_strategy) {
            // Using signed integer (int64_t) instead of size_t for the python API, since it can be useful to
            // return -1 as the virtual boundary when inspecting the graph.
            if (node1 < 0 || node2 < 0)
                throw std::invalid_argument("Node indices must be non-negative.");
            std::vector<size_t> observables_vec(observables.begin(), observables.end());
            self.add_or_merge_edge(
                node1, node2, observables_vec, weight, error_probability, merge_strategy_from_string(merge_strategy));
        },
        "node1"_a,
        "node2"_a,
        "observables"_a,
        "weight"_a,
        "error_probability"_a,
        "merge_strategy"_a);
    g.def(
        "add_boundary_edge",
        [](pm::UserGraph &self,
           int64_t node,
           const std::set<size_t> &observables,
           double weight,
           double error_probability,
           const std::string &merge_strategy) {
            // Using signed integer (int64_t) instead of size_t for the python API, since it can be useful to
            // return -1 as the virtual boundary when inspecting the graph.
            if (node < 0)
                throw std::invalid_argument("Node index must be non-negative.");
            std::vector<size_t> observables_vec(observables.begin(), observables.end());
            self.add_or_merge_boundary_edge(
                node, observables_vec, weight, error_probability, merge_strategy_from_string(merge_strategy));
        },
        "node"_a,
        "observables"_a,
        "weight"_a,
        "error_probability"_a,
        "merge_strategy"_a);
    g.def("set_boundary", &pm::UserGraph::set_boundary, "boundary"_a);
    g.def("get_boundary", &pm::UserGraph::get_boundary);
    g.def("get_num_observables", &pm::UserGraph::get_num_observables);
    g.def("set_min_num_observables", &pm::UserGraph::set_min_num_observables, "num_observables"_a);
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
    g.def("get_edges", [](const pm::UserGraph &self) {
        py::list edges;
        for (size_t i = 0; i < self.nodes.size(); i++) {
            for (size_t j = 0; j < self.nodes[i].neighbors.size(); j++) {
                auto n = self.nodes[i].neighbors[j];
                size_t v = n.node;
                if (i < v) {
                    double p;
                    if (n.error_probability < 0 || n.error_probability > 1) {
                        p = -1.0;
                    } else {
                        p = n.error_probability;
                    }
                    std::set<size_t> observables_set(n.observable_indices.begin(), n.observable_indices.end());
                    py::dict attrs("fault_ids"_a = observables_set, "weight"_a = n.weight, "error_probability"_a = p);
                    if (v == SIZE_MAX) {
                        py::tuple edge_props = py::make_tuple(i, py::none(), attrs);
                        edges.append(edge_props);
                    } else {
                        py::tuple edge_props = py::make_tuple(i, v, attrs);
                        edges.append(edge_props);
                    }
                }
            }
        }
        return edges;
    });
    g.def("has_edge", &pm::UserGraph::has_edge, "node1"_a, "node2"_a);
    g.def("has_boundary_edge", &pm::UserGraph::has_boundary_edge, "node"_a);
    g.def("get_edge_data", [](const pm::UserGraph &self, size_t node1, size_t node2) {
        if (node1 >= self.nodes.size())
            throw std::invalid_argument("node1 (" + std::to_string(node1) + ") not in graph");
        size_t idx = self.nodes[node1].index_of_neighbor(node2);
        if (idx == SIZE_MAX)
            throw std::invalid_argument(
                "Edge (" + std::to_string(node1) + ", " + std::to_string(node2) + ") not in graph.");
        auto n = self.nodes[node1].neighbors[idx];
        std::set<size_t> observables_set(n.observable_indices.begin(), n.observable_indices.end());
        py::dict attrs(
            "fault_ids"_a = observables_set, "weight"_a = n.weight, "error_probability"_a = n.error_probability);
        return attrs;
    }, "node1"_a, "node2"_a);
    g.def("get_boundary_edge_data", [](const pm::UserGraph &self, size_t node) {
            if (node >= self.nodes.size())
                throw std::invalid_argument("node (" + std::to_string(node) + ") not in graph");
            size_t idx = self.nodes[node].index_of_neighbor(SIZE_MAX);
            if (idx == SIZE_MAX)
                throw std::invalid_argument(
                    "Boundary edge (" + std::to_string(node) + ",) not in graph.");
            auto n = self.nodes[node].neighbors[idx];
            std::set<size_t> observables_set(n.observable_indices.begin(), n.observable_indices.end());
            py::dict attrs(
                "fault_ids"_a = observables_set, "weight"_a = n.weight, "error_probability"_a = n.error_probability);
            return attrs;
        }, "node"_a);
    m.def("detector_error_model_to_matching_graph", [](const char *dem_string) {
        auto dem = stim::DetectorErrorModel(dem_string);
        return pm::detector_error_model_to_user_graph(dem);
    });

    m.def(
        "sparse_column_check_matrix_to_matching_graph",
        [](const py::object &check_matrix,
           const py::array_t<double> &weights,
           const py::array_t<double> &error_probabilities,
           const std::string &merge_strategy,
           bool use_virtual_boundary_node,
           size_t num_repetitions,
           const py::array_t<double> &timelike_weights,
           const py::array_t<double> &measurement_error_probabilities) {
            py::object csc_matrix = py::module_::import("scipy.sparse").attr("csc_matrix");
            if (!py::isinstance(check_matrix, csc_matrix))
                throw std::invalid_argument("Check matrix must be a `scipy.sparse.csc_matrix`.");
            // Extract key attributes from scipy.sparse.csc_matrix:
            // https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
            // As stated in the scipy documentation for the CSC representation: "The row indices for column i
            // are stored in indices[indptr[i]:indptr[i+1]] and their corresponding values are stored in
            // data[indptr[i]:indptr[i+1]]."
            auto data_arr = check_matrix.attr("data").cast<py::array_t<uint8_t>>();
            auto data = data_arr.unchecked<1>();
            auto indices_arr = check_matrix.attr("indices").cast<py::array_t<int64_t>>();
            auto indices = indices_arr.unchecked<1>();
            auto indptr_arr = check_matrix.attr("indptr").cast<py::array_t<int64_t>>();
            auto indptr = indptr_arr.unchecked<1>();

            py::tuple shape = check_matrix.attr("shape");
            size_t num_rows = shape[0].cast<size_t>();  // The number of nodes in the matching graph
            size_t num_cols = shape[1].cast<size_t>();  // The number of edges in the matching graph

            auto weights_unchecked = weights.unchecked<1>();
            // Check weights array size is correct
            if (weights_unchecked.size() != num_cols)
                throw std::invalid_argument(
                    "The size of the `weights` array (" + std::to_string(weights_unchecked.size()) +
                    ") should match the number of columns in the check matrix (" + std::to_string(num_cols) + ")");
            auto error_probabilities_unchecked = error_probabilities.unchecked<1>();
            // Check error_probabilities array is correct
            if (error_probabilities_unchecked.size() != num_cols)
                throw std::invalid_argument(
                    "The size of the `error_probabilities` array (" +
                    std::to_string(error_probabilities_unchecked.size()) +
                    ") should match the number of columns in the check matrix (" + std::to_string(num_cols) + ")");

            // Check indptr array size is correct
            if (indptr.size() != num_cols + 1)
                throw std::invalid_argument(
                    "`check_matrix.indptr` size (" + std::to_string(indptr.size()) +
                    ") must be 1 larger than number of columns (" + std::to_string(num_cols) + ").");

            // Check data is the same size as indices
            if (data.size() != indices.size())
                throw std::invalid_argument("`check_matrix.data` must be the same size as `check_matrix.indices`");

            // Check data only contains ones
            for (py::ssize_t i = 0; i < data.size(); i++) {
                if (data(i) == 0) {
                    throw std::invalid_argument(
                        "`check_matrix.data` must only contain ones, but a zero was found. First call "
                        "`check_matrix.eliminate_zeros()` before using this method.");
                } else if (data(i) != 1) {
                    throw std::invalid_argument(
                        "`check_matrix.data` must only contain ones, but the element " + std::to_string(data(i)) +
                        " was found.");
                }
            }

            auto merge_strategy_enum = merge_strategy_from_string(merge_strategy);

            // Now construct the graph
            size_t num_detectors = num_rows * num_repetitions;
            pm::UserGraph graph(num_detectors, num_cols);
            // Each column corresponds to an edge. Iterate over the columns, adding the edges to the graph.
            // Also iterate over the number of repetitions (in case num_repetitions > 1)
            for (size_t rep = 0; rep < num_repetitions; rep++) {
                for (py::ssize_t c = 0; c < num_cols; c++) {
                    auto idx_start = indptr[c];
                    auto idx_end = indptr[c + 1];
                    auto num_dets = idx_end - idx_start;
                    if (idx_start > indices.size() - 1 && idx_start != idx_end)
                        throw std::invalid_argument(
                            "`check_matrix.indptr` elements must not exceed size of `check_matrix.indices`");
                    if (num_dets == 2) {
                        graph.add_or_merge_edge(
                            indices(idx_start) + num_rows * rep,
                            indices(idx_start + 1) + num_rows * rep,
                            {(size_t)c},
                            weights_unchecked(c),
                            error_probabilities_unchecked(c),
                            merge_strategy_enum);
                    } else if (num_dets == 1) {
                        if (use_virtual_boundary_node) {
                            graph.add_or_merge_boundary_edge(
                                indices(idx_start) + num_rows * rep,
                                {(size_t)c},
                                weights_unchecked(c),
                                error_probabilities_unchecked(c),
                                merge_strategy_enum);
                        } else {
                            graph.add_or_merge_edge(
                                indices(idx_start) + num_rows * rep,
                                num_detectors,
                                {(size_t)c},
                                weights_unchecked(c),
                                error_probabilities_unchecked(c),
                                merge_strategy_enum);
                        }
                    } else {
                        throw std::invalid_argument(
                            "`check_matrix` must contain at most two ones per column, but column " + std::to_string(c) +
                            " has " + std::to_string(num_dets) + " ones.");
                    }
                }
            }

            if (num_repetitions > 1) {
                if (timelike_weights.is(py::none()))
                    throw std::invalid_argument("must provide `timelike_weights` for repetitions > 1.");
                if (measurement_error_probabilities.is(py::none()))
                    throw std::invalid_argument("must provide `measurement_error_probabilities` for repetitions > 1.");
                auto t_weights = timelike_weights.unchecked<1>();
                if (t_weights.size() != num_rows) {
                    throw std::invalid_argument(
                        "timelike_weights has length " + std::to_string(t_weights.size()) +
                        " but its length must equal the number of columns in the check matrix (" +
                        std::to_string(num_rows) + ").");
                }
                auto meas_errs = measurement_error_probabilities.unchecked<1>();
                if (meas_errs.size() != num_rows) {
                    throw std::invalid_argument(
                        "`measurement_error_probabilities` has length " + std::to_string(meas_errs.size()) +
                        " but its length must equal the number of columns in the check matrix (" +
                        std::to_string(num_rows) + ").");
                }

                for (size_t rep = 0; rep < num_repetitions - 1; rep++) {
                    for (size_t row = 0; row < num_rows; row++) {
                        graph.add_or_merge_edge(
                            row + rep * num_rows,
                            row + (rep + 1) * num_rows,
                            {},
                            t_weights(row),
                            meas_errs(row),
                            merge_strategy_enum);
                    }
                }
            }

            // Set the boundary if not using a virtual boundary and if a boundary edge was added
            if (!use_virtual_boundary_node && graph.nodes.size() == num_detectors + 1)
                graph.set_boundary({num_detectors});

            return graph;
        },
        "check_matrix"_a,
        "weights"_a,
        "error_probabilities"_a,
        "merge_strategy"_a = "smallest-weight",
        "use_virtual_boundary_node"_a = false,
        "num_repetitions"_a = 1,
        "timelike_weights"_a = py::none(),
        "measurement_error_probabilities"_a = py::none());
}

//// Convert match edges to a vector of int64_t
//for (auto& e : mwpm.flooder.match_edges) {
//    match_edges.push_back(e.loc_from - &mwpm.flooder.graph.nodes[0]);
//    if (e.loc_to) {
//        match_edges.push_back(e.loc_to - &mwpm.flooder.graph.nodes[0]);
//    } else {
//        match_edges.push_back(-1);
//    }
//}
//
//// Put observables in a vector of uint64_t, if present
//if (num_observables <= sizeof(pm::obs_int) * 8 && return_obs_masks_if_present) {
//    for (auto& e : mwpm.flooder.match_edges) {
//        observable_masks.push_back(e.obs_mask);
//    }
//}