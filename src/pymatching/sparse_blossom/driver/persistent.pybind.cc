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

// PersistentMatcher: a matching graph built once from a check matrix (exactly as
// sparse_column_check_matrix_to_matching_graph builds it, merge strategy "smallest-weight", no virtual boundary node,
// one repetition), whose edge weights are then rewritten in place for every decode instead of
// rebuilding the graph and its flooder.
//
// Rewriting reproduces what a rebuild would produce: the graph's structure (nodes, neighbour order) does not depend on
// the weights; per user edge the kept column is the first column with the strictly smallest weight (the smallest-weight
// merge rule applied in column order); the discretisation is UserGraph::iter_discretized_edges with the normalising
// constant of UserGraph::get_edge_weight_normalising_constant; negative weights update the same sets as
// MatchingGraph::add_edge / add_boundary_edge, and the flooder's copies are refreshed with
// sync_negative_weight_observables_and_detection_events, as UserGraph::to_mwpm does.
//
// More than 64 observables (e.g. the identity faults matrix, so that decode returns the full correction as
// Matching.decode does): the Mwpm then carries a search graph, built as UserGraph::to_mwpm(..., true) builds it, whose
// weights, observable lists and sign markers are rewritten in place too; the matching graph's observable masks stay 0
// as MatchingGraph::add_edge leaves them.
#include "pymatching/sparse_blossom/driver/persistent.pybind.h"

#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "pybind11/numpy.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/user_graph.h"

namespace py = pybind11;
using namespace py::literals;

namespace {

struct Slot {
    size_t node;
    size_t index;
};

class PersistentMatcher {
   public:
    PersistentMatcher(
        const py::array_t<int64_t> &h_indptr,
        const py::array_t<int64_t> &h_indices,
        size_t num_rows,
        size_t num_cols,
        const py::array_t<int64_t> &f_indptr,
        const py::array_t<int64_t> &f_indices,
        size_t num_observables,
        const py::array_t<double> &weights)
        : num_detectors(num_rows), num_cols(num_cols), big_obs(num_observables > sizeof(pm::obs_int) * 8) {
        auto hp = h_indptr.unchecked<1>();
        auto hi = h_indices.unchecked<1>();
        auto fp = f_indptr.unchecked<1>();
        auto fi = f_indices.unchecked<1>();
        auto w = weights.unchecked<1>();
        if ((size_t)w.size() != num_cols)
            throw std::invalid_argument("weights must have one entry per column");

        pm::UserGraph graph(num_detectors, num_observables);
        col_obs.resize(num_cols);
        col_obs_mask.assign(num_cols, 0);
        col_nodes.resize(num_cols);
        for (size_t c = 0; c < num_cols; c++) {
            for (auto q = fp(c); q < fp(c + 1); q++) {
                col_obs[c].push_back((size_t)fi(q));
                if (!big_obs)
                    col_obs_mask[c] ^= (pm::obs_int)1 << fi(q);
            }
            auto n = hp(c + 1) - hp(c);
            if (n == 2) {
                col_nodes[c] = {(size_t)hi(hp(c)), (size_t)hi(hp(c) + 1)};
            } else if (n == 1) {
                col_nodes[c] = {(size_t)hi(hp(c)), num_detectors};
            } else if (n != 0) {
                throw std::invalid_argument("check matrix must contain at most two ones per column");
            } else {
                col_nodes[c] = {SIZE_MAX, SIZE_MAX};
                continue;
            }
            graph.add_or_merge_edge(
                col_nodes[c].first, col_nodes[c].second, col_obs[c], w(c), -1, pm::SMALLEST_WEIGHT);
        }
        if (graph.nodes.size() == num_detectors + 1)
            graph.set_boundary({num_detectors});
        has_boundary = graph.nodes.size() == num_detectors + 1;

        // user edges in graph order, and the columns of each (in column order)
        std::map<std::pair<size_t, size_t>, size_t> edge_of_pair;
        for (auto &e : graph.edges) {
            edge_of_pair[{e.node1, e.node2}] = edge_nodes.size();
            edge_nodes.push_back({e.node1, e.node2});
        }
        edge_cols.resize(edge_nodes.size());
        for (size_t c = 0; c < num_cols; c++) {
            if (col_nodes[c].first == SIZE_MAX)
                continue;
            auto it = edge_of_pair.find(col_nodes[c]);
            if (it == edge_of_pair.end())
                it = edge_of_pair.find({col_nodes[c].second, col_nodes[c].first});
            edge_cols[it->second].push_back(c);
        }

        mwpm = graph.to_mwpm(pm::NUM_DISTINCT_WEIGHTS, false);
        auto &g = mwpm.flooder.graph;
        auto &sg = mwpm.search_flooder.graph;
        // locate each edge's weight slots in the matching graph (and in the search graph)
        for (auto &uv : edge_nodes) {
            bool b1 = is_boundary(uv.first), b2 = is_boundary(uv.second);
            std::vector<Slot> slots, sslots;
            if (b1 || b2) {
                size_t u = b1 ? uv.second : uv.first;
                slots.push_back({u, find_neighbor(g, u, nullptr)});
                if (big_obs)
                    sslots.push_back({u, find_neighbor(sg, u, nullptr)});
            } else {
                slots.push_back({uv.first, find_neighbor(g, uv.first, &g.nodes[uv.second])});
                slots.push_back({uv.second, find_neighbor(g, uv.second, &g.nodes[uv.first])});
                if (big_obs) {
                    sslots.push_back({uv.first, find_neighbor(sg, uv.first, &sg.nodes[uv.second])});
                    sslots.push_back({uv.second, find_neighbor(sg, uv.second, &sg.nodes[uv.first])});
                }
            }
            edge_slots.push_back(slots);
            edge_sslots.push_back(sslots);
        }
        edge_weight.assign(edge_nodes.size(), 0.0);
        edge_kept.assign(edge_nodes.size(), 0);
        set_weights(weights);
    }

    void set_weights(const py::array_t<double> &weights) {
        auto w = weights.unchecked<1>();
        if ((size_t)w.size() != num_cols)
            throw std::invalid_argument("weights must have one entry per column");
        // smallest-weight merge, in column order
        double max_abs = 0;
        bool all_integral = true;
        for (size_t e = 0; e < edge_nodes.size(); e++) {
            auto &cols = edge_cols[e];
            size_t kept = cols[0];
            double we = w(kept);
            for (size_t k = 1; k < cols.size(); k++) {
                if (w(cols[k]) < we) {
                    we = w(cols[k]);
                    kept = cols[k];
                }
            }
            edge_weight[e] = we;
            edge_kept[e] = kept;
            if (std::abs(we) > max_abs)
                max_abs = std::abs(we);
            if (round(we) != we)
                all_integral = false;
        }
        if (max_abs > pm::MAX_USER_EDGE_WEIGHT)
            throw std::invalid_argument("maximum absolute edge weight exceeded");
        double nc = all_integral ? 1.0 : (double)(pm::NUM_DISTINCT_WEIGHTS - 1) / max_abs;

        auto &g = mwpm.flooder.graph;
        g.negative_weight_detection_events_set.clear();
        g.negative_weight_observables_set.clear();
        g.negative_weight_sum = 0;
        for (size_t e = 0; e < edge_nodes.size(); e++) {
            pm::signed_weight_int wi = (pm::signed_weight_int)round(edge_weight[e] * nc);
            wi *= 2;
            size_t c = edge_kept[e];
            auto &slots = edge_slots[e];
            if (wi < 0) {
                g.update_negative_weight_observables(col_obs[c]);
                for (auto &s : slots)
                    g.update_negative_weight_detection_events(s.node);
                g.negative_weight_sum += wi;
            }
            for (auto &s : slots) {
                g.nodes[s.node].neighbor_weights[s.index] = std::abs(wi);
                g.nodes[s.node].neighbor_observables[s.index] = col_obs_mask[c];
            }
        }
        g.normalising_constant = nc;
        if (big_obs) {
            auto &sg = mwpm.search_flooder.graph;
            sg.negative_weight_edges.clear();
            for (size_t e = 0; e < edge_nodes.size(); e++) {
                pm::signed_weight_int wi = 2 * (pm::signed_weight_int)round(edge_weight[e] * nc);
                size_t c = edge_kept[e];
                auto &ss = edge_sslots[e];
                if (wi < 0)
                    sg.negative_weight_edges.push_back(
                        {ss[0].node, ss.size() == 2 ? ss[1].node : SIZE_MAX});
                for (auto &s : ss) {
                    sg.nodes[s.node].neighbor_weights[s.index] = std::abs(wi);
                    sg.nodes[s.node].neighbor_observable_indices[s.index] = col_obs[c];
                    sg.nodes[s.node].neighbor_markers[s.index] = wi < 0 ? pm::WEIGHT_SIGN : 0;
                }
            }
        }
        mwpm.flooder.sync_negative_weight_observables_and_detection_events();
    }

    // same call as UserGraph.decode (enable_correlations = false); returns the observables crossed
    py::array_t<uint8_t> decode(const py::array_t<uint64_t> &detection_events) {
        std::vector<uint64_t> ev(detection_events.data(), detection_events.data() + detection_events.size());
        size_t num_obs = mwpm.flooder.graph.num_observables;
        py::array_t<uint8_t> out(num_obs);
        std::fill(out.mutable_data(), out.mutable_data() + num_obs, 0);
        pm::total_weight_int weight = 0;
        pm::decode_detection_events(mwpm, ev, out.mutable_data(), weight, false);
        return out;
    }

    size_t get_num_detectors() const {
        return num_detectors;
    }

   private:
    bool is_boundary(size_t node) const {
        return node == SIZE_MAX || (has_boundary && node == num_detectors);
    }

    template <typename Graph>
    static size_t find_neighbor(Graph &g, size_t u, const void *v) {
        auto &nb = g.nodes[u].neighbors;
        for (size_t j = 0; j < nb.size(); j++)
            if ((const void *)nb[j] == v)
                return j;
        throw std::logic_error("edge slot not found");
    }

    size_t num_detectors, num_cols;
    bool big_obs;
    bool has_boundary = false;
    pm::Mwpm mwpm;
    std::vector<std::vector<size_t>> col_obs;
    std::vector<pm::obs_int> col_obs_mask;
    std::vector<std::pair<size_t, size_t>> col_nodes;
    std::vector<std::pair<size_t, size_t>> edge_nodes;
    std::vector<std::vector<size_t>> edge_cols;
    std::vector<std::vector<Slot>> edge_slots;
    std::vector<std::vector<Slot>> edge_sslots;
    std::vector<double> edge_weight;
    std::vector<size_t> edge_kept;
};

}  // namespace

void pm_pybind::pybind_persistent_matcher(py::module &m) {
    py::class_<PersistentMatcher>(m, "PersistentMatcher")
        .def(
            py::init<
                const py::array_t<int64_t> &,
                const py::array_t<int64_t> &,
                size_t,
                size_t,
                const py::array_t<int64_t> &,
                const py::array_t<int64_t> &,
                size_t,
                const py::array_t<double> &>(),
            "h_indptr"_a,
            "h_indices"_a,
            "num_rows"_a,
            "num_cols"_a,
            "f_indptr"_a,
            "f_indices"_a,
            "num_observables"_a,
            "weights"_a)
        .def("set_weights", &PersistentMatcher::set_weights, "weights"_a)
        .def("decode", &PersistentMatcher::decode, "detection_events"_a)
        .def("get_num_detectors", &PersistentMatcher::get_num_detectors);
}
