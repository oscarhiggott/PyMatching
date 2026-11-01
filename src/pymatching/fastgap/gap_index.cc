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

#include "pymatching/fastgap/gap_index.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <queue>

namespace pm {
namespace fastgap {

uint64_t fnv1a(const void* data, size_t len, uint64_t h) {
    auto p = (const uint8_t*)data;
    for (size_t i = 0; i < len; i++) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

void GapIndex::build_graph() {
    decode_user_graph = make_decode_user_graph(spec);
    pm::Mwpm mwpm = decode_user_graph.to_mwpm(pm::NUM_DISTINCT_WEIGHTS, false);
    if (mwpm.flooder.negative_weight_sum != 0 || !mwpm.flooder.negative_weight_detection_events.empty())
        throw std::invalid_argument("fastgap v1 does not support negative edge weights.");
    graph = InteriorGraph::from_matching_graph(mwpm.flooder.graph, spec.observable);
    if (graph.num_nodes < spec.num_detectors)
        throw std::logic_error("decode graph has fewer nodes than the DEM has detectors");
    h_l = graph.side_potential(SIDE_L);
    h_r = graph.side_potential(SIDE_R);

    // The exact baseline graph must share the internal units (CLAUDE.md §2.3). It also carries
    // any detector-free L-R mechanism, as a boundary edge of the observable detector.
    pm::UserGraph baseline = make_baseline_user_graph(spec);
    pm::MatchingGraph bg = baseline.to_matching_graph(pm::NUM_DISTINCT_WEIGHTS);
    if (bg.normalising_constant != graph.normalising_constant)
        throw std::invalid_argument(
            "The exact-baseline graph and fastgap's decode graph have different normalising constants (" +
            std::to_string(bg.normalising_constant) + " vs " + std::to_string(graph.normalising_constant) +
            "), so their internal units differ. This happens when a detector-free error that flips the "
            "observable is the heaviest edge; it is not supported.");
    d_lr = DIST_INF;
    size_t obs_det = spec.baseline_obs_detector();
    if (obs_det < bg.nodes.size()) {
        auto& node = bg.nodes[obs_det];
        for (size_t j = 0; j < node.neighbors.size(); j++) {
            if (node.neighbors[j] == nullptr)
                d_lr = std::min(d_lr, (dist_int)node.neighbor_weights[j]);
        }
    }
    for (size_t x = 0; x < graph.num_nodes; x++) {
        if (!is_inf(h_l[x]) && !is_inf(h_r[x]))
            d_lr = std::min(d_lr, h_l[x] + h_r[x]);
    }
}

std::shared_ptr<GapIndex> GapIndex::create_from_dem(const std::string& dem_text, size_t observable) {
    auto idx = std::make_shared<GapIndex>();
    idx->dem_text = dem_text;
    idx->spec = canonicalize_dem(stim::DetectorErrorModel(dem_text), observable);
    idx->build_graph();
    return idx;
}

std::shared_ptr<GapIndex> GapIndex::create_from_reference_edges(
    size_t num_detectors, size_t obs_detector, std::vector<RefEdge> edges) {
    auto idx = std::make_shared<GapIndex>();
    idx->spec = spec_from_reference_edges(num_detectors, obs_detector, std::move(edges));
    idx->build_graph();
    return idx;
}

void GapIndex::set_landmarks(std::vector<LandmarkSource> sources, const std::string& config) {
    size_t n = graph.num_nodes;
    for (auto& s : sources) {
        if (s.nodes.empty())
            throw std::invalid_argument("landmark '" + s.label + "' has no nodes");
        for (auto v : s.nodes) {
            if (v >= n)
                throw std::invalid_argument("landmark '" + s.label + "' refers to a node outside the graph");
        }
    }
    landmarks = std::move(sources);
    strategy_config = config;
    uint64_t h = fnv1a(config.data(), config.size());
    for (auto& s : landmarks) {
        h = fnv1a(s.label.data(), s.label.size(), h);
        h = fnv1a(s.nodes.data(), s.nodes.size() * sizeof(uint32_t), h);
    }
    strategy_hash = h;

    // Potentials, node-major.
    size_t k = potential_stride();
    potentials.assign(n * k, DIST_INF);
    for (size_t x = 0; x < n; x++) {
        potentials[x * k + 0] = h_l[x];
        potentials[x * k + 1] = h_r[x];
    }
    for (size_t l = 0; l < landmarks.size(); l++) {
        std::vector<std::pair<uint32_t, dist_int>> init;
        for (auto v : landmarks[l].nodes)
            init.push_back({v, 0});
        auto d = graph.dijkstra(init);
        for (size_t x = 0; x < n; x++)
            potentials[x * k + 2 + l] = d[x];
    }

    derive_saturated_potentials();
    finalized = true;
}

void GapIndex::derive_saturated_potentials() {
    potentials32.resize(potentials.size());
    for (size_t i = 0; i < potentials.size(); i++)
        potentials32[i] = saturate_potential(potentials[i]);
}

dist_int GapIndex::exact_distance(uint32_t u, uint32_t v, InteriorGraph::AStarScratch& scratch) const {
    if (u == v)
        return 0;
    return graph.astar(u, v, [&](uint32_t x) { return lower_bound(x, v); }, scratch);
}

pm::Mwpm GapIndex::make_mwpm() const {
    std::lock_guard<std::mutex> lock(user_graph_mutex);
    return decode_user_graph.to_mwpm(pm::NUM_DISTINCT_WEIGHTS, false);
}

size_t GapIndex::potential_bytes() const {
    return potentials.size() * sizeof(dist_int) + potentials32.size() * sizeof(int32_t);
}

// ---------------------------------------------------------------------------------------------
// Serialisation. Little-endian raw dumps; the graph is rebuilt from the source on load and its
// fingerprint is compared, so a stale or mismatched index file is rejected.

namespace {

constexpr char MAGIC[8] = {'F', 'A', 'S', 'T', 'G', 'A', 'P', '1'};
// Version 2 dropped the local distance tables.
constexpr uint32_t FORMAT_VERSION = 2;

struct Writer {
    std::ofstream out;
    explicit Writer(const std::string& path) : out(path, std::ios::binary) {
        if (!out)
            throw std::runtime_error("cannot open '" + path + "' for writing");
    }
    template <typename T>
    void pod(const T& v) {
        out.write((const char*)&v, sizeof(T));
    }
    void str(const std::string& s) {
        pod<uint64_t>(s.size());
        out.write(s.data(), s.size());
    }
    template <typename T>
    void vec(const std::vector<T>& v) {
        pod<uint64_t>(v.size());
        out.write((const char*)v.data(), v.size() * sizeof(T));
    }
};

struct Reader {
    std::ifstream in;
    explicit Reader(const std::string& path) : in(path, std::ios::binary) {
        if (!in)
            throw std::runtime_error("cannot open '" + path + "' for reading");
    }
    void check() {
        if (!in)
            throw std::runtime_error("truncated or corrupt fastgap index file");
    }
    template <typename T>
    T pod() {
        T v;
        in.read((char*)&v, sizeof(T));
        check();
        return v;
    }
    std::string str() {
        auto n = pod<uint64_t>();
        std::string s(n, '\0');
        in.read(s.data(), n);
        check();
        return s;
    }
    template <typename T>
    std::vector<T> vec() {
        auto n = pod<uint64_t>();
        std::vector<T> v(n);
        in.read((char*)v.data(), n * sizeof(T));
        check();
        return v;
    }
};

}  // namespace

void GapIndex::save(const std::string& path) const {
    if (!finalized)
        throw std::logic_error("cannot save an index before set_landmarks");
    Writer w(path);
    w.out.write(MAGIC, sizeof(MAGIC));
    w.pod<uint32_t>(FORMAT_VERSION);
    w.pod<uint8_t>(spec.kind);
    w.pod<uint64_t>(spec.observable);
    if (spec.kind == GraphSpec::FROM_DEM) {
        w.str(dem_text);
    } else {
        w.pod<uint64_t>(spec.num_detectors);
        w.pod<uint64_t>(spec.ref_obs_detector);
        w.pod<uint64_t>(spec.ref_edges.size());
        for (auto& e : spec.ref_edges) {
            w.pod<uint64_t>(e.u);
            w.pod<uint64_t>(e.v);
            w.pod<uint8_t>(e.side_l);
            w.pod<double>(e.weight);
            w.pod<double>(e.error_probability);
        }
    }
    w.pod<uint64_t>(graph.fingerprint());
    w.pod<dist_int>(d_lr);
    w.str(strategy_config);
    w.pod<uint64_t>(strategy_hash);
    w.pod<uint64_t>(landmarks.size());
    for (auto& s : landmarks) {
        w.str(s.label);
        w.vec(s.nodes);
    }
    w.vec(potentials);
    if (!w.out)
        throw std::runtime_error("error while writing '" + path + "'");
}

std::shared_ptr<GapIndex> GapIndex::load(const std::string& path) {
    Reader r(path);
    char magic[8];
    r.in.read(magic, 8);
    r.check();
    if (std::memcmp(magic, MAGIC, 8) != 0)
        throw std::runtime_error("'" + path + "' is not a fastgap index file");
    if (r.pod<uint32_t>() != FORMAT_VERSION)
        throw std::runtime_error("unsupported fastgap index format version");
    auto kind = (GraphSpec::Kind)r.pod<uint8_t>();
    auto observable = r.pod<uint64_t>();
    std::shared_ptr<GapIndex> idx;
    if (kind == GraphSpec::FROM_DEM) {
        idx = create_from_dem(r.str(), observable);
    } else {
        auto num_detectors = r.pod<uint64_t>();
        auto obs_detector = r.pod<uint64_t>();
        std::vector<RefEdge> edges(r.pod<uint64_t>());
        for (auto& e : edges) {
            e.u = r.pod<uint64_t>();
            e.v = r.pod<uint64_t>();
            e.side_l = r.pod<uint8_t>() != 0;
            e.weight = r.pod<double>();
            e.error_probability = r.pod<double>();
        }
        idx = create_from_reference_edges(num_detectors, obs_detector, std::move(edges));
    }
    if (r.pod<uint64_t>() != idx->graph.fingerprint())
        throw std::runtime_error("fastgap index file does not match the graph rebuilt from its source");
    if (r.pod<dist_int>() != idx->d_lr)
        throw std::runtime_error("fastgap index file has a different d_LR than the rebuilt graph");
    idx->strategy_config = r.str();
    idx->strategy_hash = r.pod<uint64_t>();
    auto num_landmarks = r.pod<uint64_t>();
    for (uint64_t i = 0; i < num_landmarks; i++) {
        LandmarkSource s;
        s.label = r.str();
        s.nodes = r.vec<uint32_t>();
        idx->landmarks.push_back(std::move(s));
    }
    idx->potentials = r.vec<dist_int>();
    size_t n = idx->graph.num_nodes;
    if (idx->potentials.size() != n * idx->potential_stride())
        throw std::runtime_error("fastgap index file has inconsistent array sizes");
    for (size_t x = 0; x < n; x++) {
        if (idx->potentials[x * idx->potential_stride()] != idx->h_l[x] ||
            idx->potentials[x * idx->potential_stride() + 1] != idx->h_r[x])
            throw std::runtime_error("fastgap index file has boundary potentials that do not match the graph");
    }
    idx->derive_saturated_potentials();
    idx->finalized = true;
    return idx;
}

}  // namespace fastgap
}  // namespace pm
