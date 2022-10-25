// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "pymatching/sparse_blossom/driver/user_graph.pybind.h"

#include "pybind11/pybind11.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "stim.h"

using namespace py::literals;

pm_pybind::CompressedSparseColumnCheckMatrix::CompressedSparseColumnCheckMatrix(const py::object &matrix) {
    py::object csc_matrix = py::module_::import("scipy.sparse").attr("csc_matrix");
    if (!py::isinstance(matrix, csc_matrix))
        throw std::invalid_argument("Check matrix must be a `scipy.sparse.csc_matrix`.");
    // Extract key attributes from scipy.sparse.csc_matrix:
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
    // As stated in the scipy documentation for the CSC representation: "The row indices for column i
    // are stored in indices[indptr[i]:indptr[i+1]] and their corresponding values are stored in
    // data[indptr[i]:indptr[i+1]]."
    data = matrix.attr("data").cast<py::array_t<uint8_t>>();
    indices = matrix.attr("indices").cast<py::array_t<int64_t>>();
    indptr = matrix.attr("indptr").cast<py::array_t<int64_t>>();

    py::tuple shape = matrix.attr("shape");
    num_rows = shape[0].cast<size_t>();  // The number of nodes in the matching graph
    num_cols = shape[1].cast<size_t>();  // The number of edges in the matching graph

    auto data_unchecked = data.unchecked<1>();

    // Check indptr array size is correct
    if ((size_t)indptr.size() != num_cols + (size_t)1)
        throw std::invalid_argument(
            "`matrix.indptr` size (" + std::to_string(indptr.size()) + ") must be 1 larger than number of columns (" +
            std::to_string(num_cols) + ").");

    // Check data is the same size as indices
    if (data_unchecked.size() != indices.size())
        throw std::invalid_argument("`matrix.data` must be the same size as `matrix.indices`");

    // Check data only contains ones
    for (py::ssize_t i = 0; i < data_unchecked.size(); i++) {
        if (data_unchecked(i) == 0) {
            throw std::invalid_argument(
                "`matrix.data` must only contain ones, but a zero was found. First call "
                "`matrix.eliminate_zeros()` before using this method.");
        } else if (data_unchecked(i) != 1) {
            throw std::invalid_argument(
                "`matrix.data` must only contain ones, but the element " + std::to_string(data_unchecked(i)) +
                " was found.");
        }
    }
}

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

            if (std::abs(weight) > pm::MAX_USER_EDGE_WEIGHT) {
                auto warnings = pybind11::module::import("warnings");
                warnings.attr("warn")(
                    "Weight " + std::to_string(weight) + " of edge (" + std::to_string(node1) + ", " +
                    std::to_string(node2) + ") exceeds maximum edge weight " +
                    std::to_string(pm::MAX_USER_EDGE_WEIGHT) + " and has not been added to the matching graph.");
                return;
            }
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

            if (std::abs(weight) > pm::MAX_USER_EDGE_WEIGHT) {
                auto warnings = pybind11::module::import("warnings");
                warnings.attr("warn")(
                    "Weight " + std::to_string(weight) + " of edge (" + std::to_string(node) +
                    ", None) exceeds maximum edge weight " + std::to_string(pm::MAX_USER_EDGE_WEIGHT) +
                    " and has not been added to the matching graph.");
                return;
            }
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
    g.def(
        "decode_to_matched_detection_events_array",
        [](pm::UserGraph &self, const py::array_t<uint64_t> &detection_events) {
            auto &mwpm = self.get_mwpm();
            std::vector<uint64_t> detection_events_vec(
                detection_events.data(), detection_events.data() + detection_events.size());
            pm::decode_detection_events_to_match_edges(mwpm, detection_events_vec);
            //        auto match_edges = new std::vector<int64_t>(2 * mwpm.flooder.match_edges.size());
            py::ssize_t num_edges = (py::ssize_t)(mwpm.flooder.match_edges.size());
            py::array_t<int64_t> match_edges = py::array_t<int64_t>(num_edges * 2);
            py::buffer_info buff = match_edges.request();
            int64_t *ptr = (int64_t *)buff.ptr;

            // Convert match edges to a vector of int64_t
            for (size_t i = 0; i < mwpm.flooder.match_edges.size(); i++) {
                auto &e = mwpm.flooder.match_edges[i];
                ptr[2 * i] = e.loc_from - &mwpm.flooder.graph.nodes[0];
                if (e.loc_to) {
                    ptr[2 * i + 1] = e.loc_to - &mwpm.flooder.graph.nodes[0];
                } else {
                    ptr[2 * i + 1] = -1;
                }
            }
            // Reshape the array
            match_edges.resize({(int)num_edges, 2});
            return match_edges;
        });
    g.def(
        "decode_to_matched_detection_events_dict",
        [](pm::UserGraph &self, const py::array_t<uint64_t> &detection_events) {
            auto &mwpm = self.get_mwpm();
            std::vector<uint64_t> detection_events_vec(
                detection_events.data(), detection_events.data() + detection_events.size());
            pm::decode_detection_events_to_match_edges(mwpm, detection_events_vec);
            py::dict match_dict;
            // Convert match edges to a vector of int64_t
            for (auto &e : mwpm.flooder.match_edges) {
                int64_t from_idx = e.loc_from - &mwpm.flooder.graph.nodes[0];
                if (e.loc_to) {
                    int64_t to_idx = e.loc_to - &mwpm.flooder.graph.nodes[0];
                    match_dict[py::cast(from_idx)] = to_idx;
                    // Also add reversed key-value pair
                    match_dict[py::cast(to_idx)] = from_idx;
                } else {
                    match_dict[py::cast(from_idx)] = py::none();
                }
            }
            return match_dict;
        });
    g.def("get_edges", [](const pm::UserGraph &self) {
        py::list edges;

        for (auto &e : self.edges) {
            double p;
            if (e.error_probability < 0 || e.error_probability > 1) {
                p = -1.0;
            } else {
                p = e.error_probability;
            }
            std::set<size_t> observables_set(e.observable_indices.begin(), e.observable_indices.end());
            py::dict attrs("fault_ids"_a = observables_set, "weight"_a = e.weight, "error_probability"_a = p);
            if (e.node2 == SIZE_MAX) {
                py::tuple edge_props = py::make_tuple(e.node1, py::none(), attrs);
                edges.append(edge_props);
            } else {
                py::tuple edge_props = py::make_tuple(e.node1, e.node2, attrs);
                edges.append(edge_props);
            }
        }
        return edges;
    });
    g.def("has_edge", &pm::UserGraph::has_edge, "node1"_a, "node2"_a);
    g.def("has_boundary_edge", &pm::UserGraph::has_boundary_edge, "node"_a);
    g.def(
        "get_edge_data",
        [](const pm::UserGraph &self, size_t node1, size_t node2) {
            if (node1 >= self.nodes.size())
                throw std::invalid_argument("node1 (" + std::to_string(node1) + ") not in graph");
            size_t idx = self.nodes[node1].index_of_neighbor(node2);
            if (idx == SIZE_MAX)
                throw std::invalid_argument(
                    "Edge (" + std::to_string(node1) + ", " + std::to_string(node2) + ") not in graph.");
            auto n = self.nodes[node1].neighbors[idx];
            std::set<size_t> observables_set(
                n.edge_it->observable_indices.begin(), n.edge_it->observable_indices.end());
            py::dict attrs(
                "fault_ids"_a = observables_set,
                "weight"_a = n.edge_it->weight,
                "error_probability"_a = n.edge_it->error_probability);
            return attrs;
        },
        "node1"_a,
        "node2"_a);
    g.def(
        "get_boundary_edge_data",
        [](const pm::UserGraph &self, size_t node) {
            if (node >= self.nodes.size())
                throw std::invalid_argument("node (" + std::to_string(node) + ") not in graph");
            size_t idx = self.nodes[node].index_of_neighbor(SIZE_MAX);
            if (idx == SIZE_MAX)
                throw std::invalid_argument("Boundary edge (" + std::to_string(node) + ",) not in graph.");
            auto n = self.nodes[node].neighbors[idx];
            std::set<size_t> observables_set(
                n.edge_it->observable_indices.begin(), n.edge_it->observable_indices.end());
            py::dict attrs(
                "fault_ids"_a = observables_set,
                "weight"_a = n.edge_it->weight,
                "error_probability"_a = n.edge_it->error_probability);
            return attrs;
        },
        "node"_a);
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
           const py::array_t<double> &measurement_error_probabilities,
           py::object &faults_matrix) {
            auto H = CompressedSparseColumnCheckMatrix(check_matrix);

            if (faults_matrix.is(py::none())) {
                faults_matrix =
                    py::module_::import("scipy.sparse")
                        .attr("eye")(
                            H.num_cols, "dtype"_a = py::module_::import("numpy").attr("uint8"), "format"_a = "csc");
            }
            auto F = CompressedSparseColumnCheckMatrix(faults_matrix);

            auto weights_unchecked = weights.unchecked<1>();
            // Check weights array size is correct
            if ((size_t)weights_unchecked.size() != H.num_cols)
                throw std::invalid_argument(
                    "The size of the `weights` array (" + std::to_string(weights_unchecked.size()) +
                    ") should match the number of columns in the check matrix (" + std::to_string(H.num_cols) + ")");
            auto error_probabilities_unchecked = error_probabilities.unchecked<1>();
            // Check error_probabilities array is correct
            if ((size_t)error_probabilities_unchecked.size() != H.num_cols)
                throw std::invalid_argument(
                    "The size of the `error_probabilities` array (" +
                    std::to_string(error_probabilities_unchecked.size()) +
                    ") should match the number of columns in the check matrix (" + std::to_string(H.num_cols) + ")");

            if (H.num_cols != F.num_cols)
                throw std::invalid_argument(
                    "`faults_matrix` array with shape (" + std::to_string(F.num_rows) + ", " +
                    std::to_string(F.num_cols) +
                    ") must have the same number of columns as the check matrix, which has shape (" +
                    std::to_string(H.num_rows) + ", " + std::to_string(H.num_cols) + ").");

            auto merge_strategy_enum = merge_strategy_from_string(merge_strategy);

            auto H_indptr_unchecked = H.indptr.unchecked<1>();
            auto H_indices_unchecked = H.indices.unchecked<1>();
            auto F_indptr_unchecked = F.indptr.unchecked<1>();
            auto F_indices_unchecked = F.indices.unchecked<1>();

            // Now construct the graph
            size_t num_detectors = H.num_rows * num_repetitions;
            pm::UserGraph graph(num_detectors, F.num_rows);
            // Each column corresponds to an edge. Iterate over the columns, adding the edges to the graph.
            // Also iterate over the number of repetitions (in case num_repetitions > 1)
            for (size_t rep = 0; rep < num_repetitions; rep++) {
                for (py::ssize_t c = 0; (size_t)c < H.num_cols; c++) {
                    auto idx_start = H_indptr_unchecked[c];
                    auto idx_end = H_indptr_unchecked[c + 1];
                    auto num_dets = idx_end - idx_start;
                    if (idx_start > H_indices_unchecked.size() - 1 && idx_start != idx_end)
                        throw std::invalid_argument(
                            "`check_matrix.indptr` elements must not exceed size of `check_matrix.indices`");
                    auto f_idx_start = F_indptr_unchecked[c];
                    auto f_idx_end = F_indptr_unchecked[c + 1];
                    std::vector<size_t> obs;
                    obs.reserve(f_idx_end - f_idx_start);
                    for (auto q = f_idx_start; q < f_idx_end; q++)
                        obs.push_back((size_t)F_indices_unchecked(q));
                    if (num_dets == 2) {
                        graph.add_or_merge_edge(
                            H_indices_unchecked(idx_start) + H.num_rows * rep,
                            H_indices_unchecked(idx_start + 1) + H.num_rows * rep,
                            obs,
                            weights_unchecked(c),
                            error_probabilities_unchecked(c),
                            merge_strategy_enum);
                    } else if (num_dets == 1) {
                        if (use_virtual_boundary_node) {
                            graph.add_or_merge_boundary_edge(
                                H_indices_unchecked(idx_start) + H.num_rows * rep,
                                obs,
                                weights_unchecked(c),
                                error_probabilities_unchecked(c),
                                merge_strategy_enum);
                        } else {
                            graph.add_or_merge_edge(
                                H_indices_unchecked(idx_start) + H.num_rows * rep,
                                num_detectors,
                                obs,
                                weights_unchecked(c),
                                error_probabilities_unchecked(c),
                                merge_strategy_enum);
                        }
                    } else if (num_dets != 0) {
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
                if ((size_t)t_weights.size() != H.num_rows) {
                    throw std::invalid_argument(
                        "timelike_weights has length " + std::to_string(t_weights.size()) +
                        " but its length must equal the number of columns in the check matrix (" +
                        std::to_string(H.num_rows) + ").");
                }
                auto meas_errs = measurement_error_probabilities.unchecked<1>();
                if ((size_t)meas_errs.size() != H.num_rows) {
                    throw std::invalid_argument(
                        "`measurement_error_probabilities` has length " + std::to_string(meas_errs.size()) +
                        " but its length must equal the number of columns in the check matrix (" +
                        std::to_string(H.num_rows) + ").");
                }

                for (size_t rep = 0; rep < num_repetitions - 1; rep++) {
                    for (size_t row = 0; row < H.num_rows; row++) {
                        graph.add_or_merge_edge(
                            row + rep * H.num_rows,
                            row + (rep + 1) * H.num_rows,
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
        "measurement_error_probabilities"_a = py::none(),
        "faults_matrix"_a = py::none());
}