// Copyright 2026 fastgap contributors
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

// Low-level bindings for the `fastgap` Python package. Everything here works in PyMatching's
// internal integer units; the Python layer converts to user units at the API boundary.

#include <set>

#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pymatching/fastgap/exact.h"
#include "pymatching/fastgap/gap_decoder.h"
#include "pymatching/fastgap/gap_index.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"

namespace py = pybind11;
using namespace py::literals;
using namespace pm::fastgap;

namespace {

template <typename T>
py::array_t<T> to_array(const std::vector<T>& v) {
    py::array_t<T> out(v.size());
    std::copy(v.begin(), v.end(), out.mutable_data());
    return out;
}

std::vector<uint64_t> to_events(const py::array_t<uint64_t, py::array::c_style | py::array::forcecast>& a) {
    return std::vector<uint64_t>(a.data(), a.data() + a.size());
}

Pricing parse_pricing(const std::string& s) {
    if (s == "index")
        return Pricing::INDEX;
    if (s == "exact")
        return Pricing::EXACT;
    throw std::invalid_argument("pricing must be 'index' or 'exact'");
}

DecodeOptions make_options(
    int64_t threshold, bool upper_bound, const std::string& pricing, bool check, bool keep_walk) {
    DecodeOptions o;
    o.threshold = threshold;
    o.upper_bound = upper_bound;
    o.pricing = parse_pricing(pricing);
    o.check = check;
    o.keep_walk = keep_walk;
    return o;
}

/// Validates a 2D uint8 bit-packed shot array and returns (num_shots, bytes_per_shot).
std::pair<size_t, size_t> check_shots(const py::array_t<uint8_t, py::array::c_style>& shots, size_t num_detectors) {
    if (shots.ndim() != 2)
        throw std::invalid_argument("shots must be a 2D uint8 array of bit-packed detection events");
    size_t bytes = shots.shape(1);
    if (bytes * 8 < num_detectors)
        throw std::invalid_argument(
            "bit-packed shots need at least " + std::to_string((num_detectors + 7) / 8) + " bytes per row");
    return {(size_t)shots.shape(0), bytes};
}

py::dict gap_results_to_dict(const std::vector<GapResult>& rs) {
    size_t n = rs.size();
    py::array_t<uint64_t> prediction(n);
    py::array_t<int64_t> w_star(n), gap_lb(n), gap_ub(n), hop_slack(n), t_decode(n), t_gap(n);
    py::array_t<bool> censored(n), simple(n), exact(n), bpass(n);
    py::array_t<uint32_t> nd(n), nb(n), nh(n);
    for (size_t i = 0; i < n; i++) {
        auto& r = rs[i];
        prediction.mutable_at(i) = r.prediction;
        w_star.mutable_at(i) = r.w_star;
        gap_lb.mutable_at(i) = is_inf(r.gap_lb) ? -1 : r.gap_lb;
        gap_ub.mutable_at(i) = is_inf(r.gap_ub) ? -1 : r.gap_ub;
        hop_slack.mutable_at(i) = r.hop_slack;
        t_decode.mutable_at(i) = r.t_decode_ns;
        t_gap.mutable_at(i) = r.t_gap_ns;
        censored.mutable_at(i) = r.censored;
        simple.mutable_at(i) = r.walk_simple;
        exact.mutable_at(i) = r.all_hops_exact;
        bpass.mutable_at(i) = r.blossom_pass;
        nd.mutable_at(i) = r.num_defects;
        nb.mutable_at(i) = r.num_blossoms;
        nh.mutable_at(i) = r.num_hops;
    }
    return py::dict(
        "prediction_mask"_a = prediction,
        "w_star"_a = w_star,
        "gap_lb"_a = gap_lb,
        "gap_ub"_a = gap_ub,
        "hop_slack"_a = hop_slack,
        "t_decode_ns"_a = t_decode,
        "t_gap_ns"_a = t_gap,
        "censored"_a = censored,
        "walk_simple"_a = simple,
        "all_hops_exact"_a = exact,
        "blossom_pass"_a = bpass,
        "num_defects"_a = nd,
        "num_blossoms"_a = nb,
        "num_hops"_a = nh);
}

py::dict exact_results_to_dict(const std::vector<ExactGap>& rs) {
    size_t n = rs.size();
    py::array_t<uint8_t> prediction(n);
    py::array_t<int64_t> w_star(n), gap(n), w_even(n), w_odd(n), t(n), t_mstar(n), t_comp(n);
    for (size_t i = 0; i < n; i++) {
        prediction.mutable_at(i) = rs[i].prediction;
        w_star.mutable_at(i) = rs[i].w_star;
        gap.mutable_at(i) = is_inf(rs[i].gap) ? -1 : rs[i].gap;
        w_even.mutable_at(i) = is_inf(rs[i].w_even) ? -1 : rs[i].w_even;
        w_odd.mutable_at(i) = is_inf(rs[i].w_odd) ? -1 : rs[i].w_odd;
        t.mutable_at(i) = rs[i].t_exact_ns;
        t_mstar.mutable_at(i) = rs[i].t_mstar_ns;
        t_comp.mutable_at(i) = rs[i].t_comp_ns;
    }
    return py::dict(
        "prediction"_a = prediction,
        "w_star"_a = w_star,
        "gap"_a = gap,
        "w_even"_a = w_even,
        "w_odd"_a = w_odd,
        "t_exact_ns"_a = t,
        "t_mstar_ns"_a = t_mstar,
        "t_comp_ns"_a = t_comp);
}

/// Compares fastgap's decode graph with the graph PyMatching builds from a reference gap DEM
/// (CLAUDE.md §7.3). Returns an empty string when edges, weights and units are identical.
std::string compare_with_reference_dem(const GapIndex& idx, const std::string& gap_dem_text) {
    if (idx.spec.kind != GraphSpec::FROM_REFERENCE_EDGES)
        return "index was not built from reference edges";
    size_t obs_det = idx.spec.ref_obs_detector;
    pm::Mwpm ref = pm::detector_error_model_to_mwpm(stim::DetectorErrorModel(gap_dem_text), pm::NUM_DISTINCT_WEIGHTS);
    auto& rg = ref.flooder.graph;
    auto& g = idx.graph;
    if (rg.normalising_constant != g.normalising_constant)
        return "normalising constants differ";
    if (rg.nodes.size() != g.num_nodes)
        return "node counts differ";
    const pm::DetectorNode* base = rg.nodes.data();
    for (size_t u = 0; u < g.num_nodes; u++) {
        // (neighbour, weight, side) multiset; neighbour -1 = boundary.
        std::multiset<std::tuple<int64_t, int64_t, int>> a, b;
        auto& node = rg.nodes[u];
        for (size_t j = 0; j < node.neighbors.size(); j++) {
            int64_t w = node.neighbor_weights[j];
            if (node.neighbors[j] == nullptr) {
                a.insert({-1, w, u == obs_det ? 0 : (int)SIDE_R});
            } else {
                size_t v = node.neighbors[j] - base;
                if (u == obs_det)
                    continue;  // its edges are compared from the other endpoint
                if (v == obs_det)
                    a.insert({-1, w, (int)SIDE_L});
                else
                    a.insert({(int64_t)v, w, 0});
            }
        }
        for (uint32_t e = g.offsets[u]; e < g.offsets[u + 1]; e++)
            b.insert({(int64_t)g.neighbors[e], g.weights[e], 0});
        if (g.boundary_side[u] != SIDE_NONE)
            b.insert({-1, g.boundary_weight[u], u == obs_det ? 0 : (int)g.boundary_side[u]});
        if (a != b)
            return "edges at detector D" + std::to_string(u) + " differ";
    }
    return "";
}

}  // namespace

PYBIND11_MODULE(_cpp_fastgap, m) {
    m.doc() = "Low-level C++ core of fastgap (internal integer units).";

    py::register_exception<BoundarylessLogicalError>(m, "BoundarylessLogicalError", PyExc_ValueError);
    py::register_exception<ParallelEdgeConflictError>(m, "ParallelEdgeConflictError", PyExc_ValueError);
    py::register_exception<InvariantViolation>(m, "InvariantViolation", PyExc_AssertionError);

    m.attr("DIST_INF") = DIST_INF;

    py::class_<GapIndex, std::shared_ptr<GapIndex>>(m, "GapIndexCore")
        .def_static("create_from_dem", &GapIndex::create_from_dem, "dem_text"_a, "observable"_a)
        .def_static(
            "create_from_reference_edges",
            [](size_t num_detectors, size_t obs_detector, const std::vector<std::tuple<int64_t, int64_t, bool, double, double>>& edges) {
                std::vector<RefEdge> out;
                for (auto& [u, v, side_l, w, p] : edges)
                    out.push_back({u < 0 ? SIZE_MAX : (size_t)u, v < 0 ? SIZE_MAX : (size_t)v, side_l, w, p});
                return GapIndex::create_from_reference_edges(num_detectors, obs_detector, std::move(out));
            },
            "num_detectors"_a,
            "obs_detector"_a,
            "edges"_a)
        .def_static("load", &GapIndex::load, "path"_a)
        .def("save", &GapIndex::save, "path"_a)
        .def(
            "set_landmarks",
            [](GapIndex& self, const std::vector<std::pair<std::string, std::vector<uint32_t>>>& sources,
               const std::string& config) {
                std::vector<LandmarkSource> ls;
                for (auto& [label, nodes] : sources)
                    ls.push_back({label, nodes});
                self.set_landmarks(std::move(ls), config);
            },
            "sources"_a,
            "config"_a)
        .def_property_readonly("kind", [](const GapIndex& s) { return s.spec.kind == GraphSpec::FROM_DEM ? "dem" : "reference"; })
        .def_property_readonly("dem_text", [](const GapIndex& s) { return s.dem_text; })
        .def_property_readonly("num_detectors", [](const GapIndex& s) { return s.spec.num_detectors; })
        .def_property_readonly("num_nodes", [](const GapIndex& s) { return s.graph.num_nodes; })
        .def_property_readonly("num_observables", [](const GapIndex& s) { return s.spec.num_observables; })
        .def_property_readonly("observable", [](const GapIndex& s) { return s.spec.observable; })
        .def_property_readonly("ref_obs_detector", [](const GapIndex& s) -> int64_t {
            return s.spec.ref_obs_detector == SIZE_MAX ? -1 : (int64_t)s.spec.ref_obs_detector;
        })
        .def_property_readonly("num_hyperedges_ignored", [](const GapIndex& s) { return s.spec.num_hyperedges_ignored; })
        .def_property_readonly("normalising_constant", [](const GapIndex& s) { return s.graph.normalising_constant; })
        .def_property_readonly("d_lr", [](const GapIndex& s) { return s.d_lr; })
        .def_property_readonly("finalized", [](const GapIndex& s) { return s.finalized; })
        .def_property_readonly("strategy_config", [](const GapIndex& s) { return s.strategy_config; })
        .def_property_readonly("strategy_hash", [](const GapIndex& s) { return s.strategy_hash; })
        .def_property_readonly("fingerprint", [](const GapIndex& s) { return s.graph.fingerprint(); })
        .def_property_readonly("potential_bytes", &GapIndex::potential_bytes)
        .def_property_readonly("median_edge_weight", [](const GapIndex& s) { return s.graph.median_edge_weight(); })
        .def_property_readonly("gauge", [](const GapIndex& s) { return to_array(s.spec.gauge); })
        .def_property_readonly("h_l", [](const GapIndex& s) { return to_array(s.h_l); })
        .def_property_readonly("h_r", [](const GapIndex& s) { return to_array(s.h_r); })
        .def_property_readonly(
            "landmarks",
            [](const GapIndex& s) {
                py::list out;
                for (auto& l : s.landmarks)
                    out.append(py::make_tuple(l.label, to_array(l.nodes)));
                return out;
            })
        .def_property_readonly(
            "potentials",
            [](const GapIndex& s) {
                size_t k = s.potential_stride();
                py::array_t<int64_t> out({s.graph.num_nodes, k});
                std::copy(s.potentials.begin(), s.potentials.end(), out.mutable_data());
                return out;
            })
        .def(
            "graph_arrays",
            [](const GapIndex& s) {
                return py::dict(
                    "offsets"_a = to_array(s.graph.offsets),
                    "neighbors"_a = to_array(s.graph.neighbors),
                    "weights"_a = to_array(s.graph.weights),
                    "boundary_weight"_a = to_array(s.graph.boundary_weight),
                    "boundary_side"_a = to_array(s.graph.boundary_side));
            })
        .def("canonical_dem_text", [](const GapIndex& s) {
            if (s.spec.kind != GraphSpec::FROM_DEM)
                throw std::invalid_argument("index was built from reference edges, not a DEM");
            return s.spec.canonical_dem.str();
        })
        .def("baseline_dem_text", [](const GapIndex& s) { return make_baseline_dem(s.spec).str(); })
        .def(
            "dijkstra",
            [](const GapIndex& s, const std::vector<std::pair<uint32_t, int64_t>>& init, int64_t cap) {
                std::vector<dist_int> d;
                {
                    py::gil_scoped_release release;
                    d = s.graph.dijkstra(init, cap < 0 ? DIST_INF : cap);
                }
                return to_array(d);
            },
            "init"_a,
            "cap"_a = -1)
        .def(
            "exact_distance",
            [](const GapIndex& s, uint32_t u, uint32_t v) {
                InteriorGraph::AStarScratch scratch;
                return s.exact_distance(u, v, scratch);
            },
            "u"_a,
            "v"_a)
        .def("lower_bound", &GapIndex::lower_bound, "u"_a, "v"_a)
        .def(
            "d_hat",
            [](const GapIndex& s, uint32_t u, uint32_t v) {
                bool exact;
                dist_int d = s.d_hat(u, v, exact);
                return py::make_tuple(d, exact);
            },
            "u"_a,
            "v"_a)
        .def("compare_with_reference_dem", &compare_with_reference_dem, "gap_dem_text"_a)
        .def(
            "dual_state",
            [](const GapIndex& s, const py::array_t<uint64_t, py::array::c_style | py::array::forcecast>& events) {
                // The decoder's final dual state for one shot (for inspection and figures).
                pm::Mwpm mwpm = s.make_mwpm();
                DualState ds;
                DualStateScratch scratch;
                extract_dual_state(mwpm, to_events(events), s.spec.observable, ds, scratch);
                py::list ancestors;
                for (size_t i = 0; i < ds.num_defects(); i++) {
                    py::list a;
                    for (uint32_t k = ds.anc_offset[i]; k < ds.anc_offset[i + 1]; k++)
                        a.append(ds.anc_blossom[k]);
                    ancestors.append(a);
                }
                return py::dict(
                    "defects"_a = to_array(ds.defects),
                    "y_leaf"_a = to_array(ds.y_leaf),
                    "rho"_a = to_array(ds.rho),
                    "y_blossom"_a = to_array(ds.y_blossom),
                    "ancestors"_a = ancestors,
                    "mate"_a = to_array(ds.mate),
                    "w_star"_a = ds.w_star,
                    "mate_l"_a = MATE_L,
                    "mate_r"_a = MATE_R);
            },
            "detection_events"_a);

    py::class_<GapDecoder>(m, "GapDecoderCore")
        .def(
            py::init([](std::shared_ptr<GapIndex> idx, size_t num_threads, bool pin) {
                return new GapDecoder(idx, num_threads, pin);
            }),
            "index"_a,
            "num_threads"_a = 1,
            "pin_threads"_a = false)
        .def_property_readonly("num_threads", &GapDecoder::num_threads)
        .def_readwrite("parallel_grain", &GapDecoder::parallel_grain)
        .def(
            "decode",
            [](GapDecoder& self, const py::array_t<uint64_t, py::array::c_style | py::array::forcecast>& events,
               int64_t threshold, bool upper_bound, const std::string& pricing, bool check, bool keep_walk) {
                auto ev = to_events(events);
                auto o = make_options(threshold, upper_bound, pricing, check, keep_walk);
                GapResult r;
                {
                    py::gil_scoped_release release;
                    r = self.decode(ev, o);
                }
                py::dict d = gap_results_to_dict({r});
                py::dict out;
                for (auto item : d)
                    out[item.first] = item.second.cast<py::array>()[py::int_(0)];
                py::list walk;
                for (auto& h : r.walk)
                    walk.append(py::make_tuple(h.from, h.to, h.price, h.exact));
                out["walk"] = walk;
                return out;
            },
            "detection_events"_a,
            "threshold"_a = -1,
            "upper_bound"_a = false,
            "pricing"_a = "index",
            "check"_a = CHECKS_DEFAULT_ON,
            "keep_walk"_a = false)
        .def(
            "decode_batch",
            [](GapDecoder& self, const py::array_t<uint8_t, py::array::c_style>& shots, int64_t threshold,
               bool upper_bound, const std::string& pricing, bool check) {
                auto [num_shots, bytes] = check_shots(shots, self.index().num_detectors());
                auto o = make_options(threshold, upper_bound, pricing, check, false);
                std::vector<GapResult> rs;
                {
                    py::gil_scoped_release release;
                    self.decode_batch(shots.data(), num_shots, bytes, o, rs);
                }
                return gap_results_to_dict(rs);
            },
            "shots"_a,
            "threshold"_a = -1,
            "upper_bound"_a = false,
            "pricing"_a = "index",
            "check"_a = CHECKS_DEFAULT_ON);

    py::class_<ExactBaseline>(m, "ExactBaselineCore")
        .def(
            py::init([](std::shared_ptr<GapIndex> idx, size_t num_threads) { return new ExactBaseline(idx, num_threads); }),
            "index"_a,
            "num_threads"_a = 1)
        .def_property_readonly("normalising_constant", &ExactBaseline::normalising_constant)
        .def(
            "decode",
            [](ExactBaseline& self, const py::array_t<uint64_t, py::array::c_style | py::array::forcecast>& events) {
                auto r = self.decode(to_events(events));
                py::dict d = exact_results_to_dict({r});
                py::dict out;
                for (auto item : d)
                    out[item.first] = item.second.cast<py::array>()[py::int_(0)];
                return out;
            },
            "detection_events"_a)
        .def(
            "decode_batch",
            [](ExactBaseline& self, const py::array_t<uint8_t, py::array::c_style>& shots, size_t num_detectors) {
                auto [num_shots, bytes] = check_shots(shots, num_detectors);
                std::vector<ExactGap> rs;
                {
                    py::gil_scoped_release release;
                    self.decode_batch(shots.data(), num_shots, bytes, rs);
                }
                return exact_results_to_dict(rs);
            },
            "shots"_a,
            "num_detectors"_a);

    m.def(
        "brute_force_gap",
        [](const GapIndex& idx, const py::array_t<uint64_t, py::array::c_style | py::array::forcecast>& events) {
            auto r = brute_force_gap(idx, to_events(events));
            return py::dict(
                "prediction"_a = r.prediction,
                "w_star"_a = r.w_star,
                "gap"_a = is_inf(r.gap) ? -1 : r.gap,
                "w_even"_a = is_inf(r.w_even) ? -1 : r.w_even,
                "w_odd"_a = is_inf(r.w_odd) ? -1 : r.w_odd);
        },
        "index"_a,
        "detection_events"_a);

    m.def(
        "pin_current_thread",
        [](size_t cpu) { return ThreadPool::pin_current_thread(cpu); },
        "cpu"_a);
}
