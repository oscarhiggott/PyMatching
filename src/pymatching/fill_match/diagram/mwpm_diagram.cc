#include "pymatching/fill_match/diagram/mwpm_diagram.h"
#include "pymatching/fill_match/driver/mwpm_decoding.h"

using namespace pm;

std::vector<std::pair<float, float>> pm::dem_detector_coords(const stim::DetectorErrorModel &dem) {
    size_t n = dem.count_detectors();
    std::set<uint64_t> all_dets;
    for (uint64_t k = 0; k < n; k++) {
        all_dets.insert(k);
    }
    auto coords = dem.get_detector_coordinates(all_dets);
    std::vector<std::pair<float, float>> result;
    for (uint64_t k = 0; k < n; k++) {
        auto p = coords.find(k);
        if (p == coords.end() || p->second.empty()) {
            result.push_back({k, 0});
            continue;
        }
        auto &cs = p->second;
        if (cs.size() == 1) {
            result.push_back({cs[0], 0});
        } else {
            result.push_back({cs[0], cs[1]});
        }

        // Do an arbitrary orthographic projection of the other axes.
        double s = 1;
        for (size_t d = 2; d < cs.size(); d++) {
            s *= 0.66;
            result.back().first += cs[d] * s;
            result.back().second += cs[d] * s / (d + 1);
        }
    }
    return result;
}

struct StateHelper {
    const Mwpm &mwpm;
    const std::vector<DetectorNode> &ns;
    const cumulative_time_int &t;
    std::vector<std::pair<float, float>> coords;
    float scale;
    std::ostream &out;

    std::pair<float, float> neighbor_coords(const DetectorNode &src, const DetectorNode *dst) {
        size_t k = &src - &ns[0];
        if (dst != nullptr) {
            size_t k2 = dst - &ns[0];
            return coords[k2];
        } else {
            return {coords[k].first + 1.0f, coords[k].second + 0.7f};
        }
    };

    std::pair<float, float> lerp_pos(const DetectorNode &src, size_t local_neighbor_index, float ratio) {
        if (ratio > 1) {
            ratio = 1;
        }
        auto c = coords[&src - &ns[0]];
        auto nc = neighbor_coords(src, src.neighbors[local_neighbor_index]);
        return {nc.first * ratio + c.first * (1 - ratio),
                nc.second * ratio + c.second * (1 - ratio)};
    }

    bool is_node_in_disjoint_area_of_ancestor_of(const DetectorNode &n, GraphFillRegion *region) {
        while (region != nullptr) {
            if (n.region_that_arrived == region) {
                return true;
            }
            region = region->blossom_parent;
        }
        return false;
    }

    bool are_nodes_in_nested_regions(const DetectorNode &src, const DetectorNode &dst) {
        if (is_node_in_disjoint_area_of_ancestor_of(src, dst.region_that_arrived)) {
            return true;
        }
        if (is_node_in_disjoint_area_of_ancestor_of(dst, src.region_that_arrived)) {
            return true;
        }
        return false;
    }

    std::vector<std::vector<std::pair<float, float>>> approximate_region_polygons(GraphFillRegion *region) {
        std::map<DetectorNode *, std::set<std::pair<float, float>>> points_by_source;
        region->do_op_for_each_node_in_total_area([&](const DetectorNode *n) {
            auto &perimeter = points_by_source[n->reached_from_source];
            auto r = n->compute_wrapped_radius_within_layer_at_time(region, t);
            if (r == 0) {
                perimeter.insert(coords[n - &ns[0]]);
                return;
            }

            for (size_t nk = 0; nk < n->neighbor_weights.size(); nk++) {
               auto w = n->neighbor_weights[nk];
               float ratio = r / std::max(1.0f, (float)w);
               auto ne = n->neighbors[nk];
               if (ne == nullptr) {
                   perimeter.insert(lerp_pos(*n, nk, ratio));
                   continue;
               }
               if (ne->reached_from_source == n->reached_from_source) {
                   continue;
               }
               if (!are_nodes_in_nested_regions(*n, *ne)) {
                   auto re = ne->compute_wrapped_radius_within_layer_at_time(region, t);
                   if (r + re > w) {
                       size_t d = r + re - w;
                       if (r < d / 2) {
                           ratio = 0;
                           std::cerr << "NO\n";
                       } else {
                           ratio = (r - d / 2) / std::max(1.0f, (float) w);
                       }
                   }
               }
//               if (ratio > 1) {
//                   throw std::invalid_argument("NO");
//               }
               perimeter.insert(lerp_pos(*n, nk, ratio));
            }
        });

        std::vector<std::vector<std::pair<float, float>>> result;
        for (const auto &b : points_by_source) {
            result.push_back({});
            auto &p = result.back();
            p.insert(p.end(), b.second.begin(), b.second.end());
            auto center = coords[b.first - &ns[0]];
            std::sort(p.begin(), p.end(), [&](std::pair<float, float> a, std::pair<float, float> b) -> bool {
                return atan2f(a.second - center.second, a.first - center.first)
                     < atan2f(b.second - center.second, b.first - center.first);
            });
        }
        return result;
    }

    void fill_region(GraphFillRegion *region, const char *fill, const char *stroke) {
        for (const auto &poly : approximate_region_polygons(region)) {
            out << " <path d=\"";
            for (size_t k = 0; k < poly.size(); k++) {
                out << (k == 0 ? 'M' : 'L');
                out << poly[k].first * scale << " " << poly[k].second * scale << " ";
            }
            out << " z\" stroke=\"" << stroke << "\" fill=\"" << fill << "\"/>\n";
        }
    }

    void draw_detector_graph_edges() {
        for (size_t k = 0; k < ns.size(); k++) {
            const auto &n = ns[k];
            for (const auto *n2_ptr : n.neighbors) {
                auto nc = neighbor_coords(n, n2_ptr);
                if (n2_ptr != nullptr) {
                    size_t k2 = n2_ptr - &ns[0];
                    if (k2 > k) {
                        out << " <line x1=\"" << coords[k].first * scale
                            << "\" x2=\"" << nc.first * scale
                            << "\" y1=\"" << coords[k].second * scale
                            << "\" y2=\"" << nc.second * scale
                            << "\" stroke=\"#CCC\"/>\n";
                    }
                } else {
                    out << " <line x1=\"" << coords[k].first * scale
                        << "\" x2=\"" << nc.first * scale
                        << "\" y1=\"" << coords[k].second * scale
                        << "\" y2=\"" << nc.second * scale
                        << "\" stroke-dasharray=\"4 2"
                        << "\" stroke=\"#FAA\"/>\n";
                }
            }
        }
    }
    void draw_unexcited_detector_nodes() {
        for (size_t k = 0; k < ns.size(); k++) {
            const auto &n = ns[k];
            if (n.reached_from_source != &n) {
                out << " <circle cx=\"" << coords[k].first * scale
                    << "\" cy=\"" << coords[k].second * scale
                    << "\" r=\"" << 3
                    << "\" stroke=\"none"
                    << "\" fill=\"#CCC\"/>\n";
            }
        }
    }

    std::set<GraphFillRegion *> find_all_regions() {
        std::set<GraphFillRegion *> regions;
        for (size_t k = 0; k < ns.size(); k++) {
            const auto &n = ns[k];
            GraphFillRegion *r = n.region_that_arrived;
            while (r != nullptr) {
                regions.insert(r);
                r = r->blossom_parent;
            }
        }
        return regions;
    }

    std::set<pm::AltTreeNode *> find_all_alt_tree_nodes() {
        std::set<pm::AltTreeNode *> tree_nodes;
        for (auto *r : find_all_regions()) {
            if (r->alt_tree_node != nullptr) {
                tree_nodes.insert(r->alt_tree_node);
            }
        }
        return tree_nodes;
    }

    void draw_detection_events() {
        for (size_t k = 0; k < ns.size(); k++) {
            const auto &n = ns[k];
            if (n.reached_from_source == &n) {
                out << " <circle cx=\"" << coords[k].first * scale
                    << "\" cy=\"" << coords[k].second * scale
                    << "\" r=\"" << 3
                    << "\" stroke=\"none"
                    << "\" fill=\"red\"/>\n";
            }
        }
    }

    void draw_pending_collisions() {
        for (size_t k = 0; k < ns.size(); k++) {
            const auto &n = ns[k];
            if (n.region_that_arrived_top == nullptr) {
                continue;
            }
            auto ev = mwpm.flooder.find_next_event_at_node_returning_neighbor_index_and_time(n);
            if (ev.first != SIZE_MAX && ev.second == t) {
                auto r1 = n.local_radius();
                auto m = n.neighbors[ev.first];
                auto r2 = m == nullptr ? Varying32{0} : m->local_radius();
                if (!r1.colliding_with(r2)) {
                    continue;
                }
                if (r1.get_distance_at_time(t) + r2.get_distance_at_time(t) != n.neighbor_weights[ev.first]) {
                    continue;
                }
                if (m != nullptr && (m->region_that_arrived_top == nullptr || m->region_that_arrived_top == n.region_that_arrived_top)) {
                    continue;
                }
                auto col = lerp_pos(n, ev.first, r1.get_distance_at_time(t) / std::max(1.0f, (float)n.neighbor_weights[ev.first]));
                out << " <circle cx=\"" << col.first * scale
                    << "\" cy=\"" << col.second * scale
                    << "\" r=\"" << 3
                    << "\" stroke=\"black"
                    << "\" fill=\"orange\"/>\n";
            }
        }
    }

    void draw_match_edges() {
        for (auto *r : find_all_regions()) {
            if (r->blossom_parent_top == r && r->radius.is_frozen()) {
                size_t k1 = r->match.edge.loc_from - &ns[0];
                if (r->match.edge.loc_to == nullptr) {
                    out << " <circle cx=\"" << coords[k1].first * scale
                        << "\" cy=\"" << coords[k1].second * scale
                        << "\" r=\"" << 8
                        << "\" stroke=\"black"
                        << "\" fill=\"#FFF\"/>\n";
                } else {
                    size_t k2 = r->match.edge.loc_to - &ns[0];
                    out << " <line x1=\"" << coords[k1].first * scale
                        << "\" x2=\"" << coords[k2].first * scale
                        << "\" y1=\"" << coords[k1].second * scale
                        << "\" y2=\"" << coords[k2].second * scale
                        << "\" stroke-width=\"" << 4
                        << "\" stroke=\"#000\"/>\n";
                }
            }
        }
    }

    void draw_alt_tree_edges() {
        for (auto *tn : find_all_alt_tree_nodes()) {
            if (tn->inner_region != nullptr) {
                size_t k1 = tn->inner_to_outer_edge.loc_from - &ns[0];
                size_t k2 = tn->inner_to_outer_edge.loc_to - &ns[0];
                out << " <line x1=\"" << coords[k1].first * scale
                    << "\" x2=\"" << coords[k2].first * scale
                    << "\" y1=\"" << coords[k1].second * scale
                    << "\" y2=\"" << coords[k2].second * scale
                    << "\" stroke=\"#000\"/>\n";
            }
            for (size_t k = 0; k < tn->children.size(); k++) {
                size_t k1 = tn->children[k].edge.loc_from - &ns[0];
                size_t k2 = tn->children[k].edge.loc_to - &ns[0];
                out << " <line x1=\"" << coords[k1].first * scale
                    << "\" x2=\"" << coords[k2].first * scale
                    << "\" y1=\"" << coords[k1].second * scale
                    << "\" y2=\"" << coords[k2].second * scale
                    << "\" stroke-dasharray=\"4 2"
                    << "\" stroke=\"#000\"/>\n";
            }
        }
    }

    void draw_blossom_edges() {
        for (auto *r : find_all_regions()) {
            for (const auto &e : r->blossom_children) {
                size_t k1 = e.edge.loc_from - &ns[0];
                size_t k2 = e.edge.loc_to - &ns[0];
                out << " <line x1=\"" << coords[k1].first * scale
                    << "\" x2=\"" << coords[k2].first * scale
                    << "\" y1=\"" << coords[k1].second * scale
                    << "\" y2=\"" << coords[k2].second * scale
                    << "\" stroke-width=\"" << 4
                    << "\" stroke=\"#22F\"/>\n";
            }
        }
    }
};

void pm::write_mwpm_svg_diagram(const std::vector<std::pair<float, float>> coords, const pm::Mwpm &mwpm, MwpmEvent focused_event, std::ostream &out) {
    out << R"SVG(<svg width="2000" height="800" version="1.1" xmlns="http://www.w3.org/2000/svg">)SVG" << '\n';

    out << " <path d=\"M 0 0 L 2000 0 L 2000 800 L 0 800 Z\" fill=\"white\" stroke=\"#CCC\"/>\n";

    StateHelper s{
        mwpm,
        mwpm.flooder.graph.nodes,
        mwpm.flooder.queue.cur_time,
        coords,
        20,
        out,
    };

    s.draw_detector_graph_edges();
    s.draw_unexcited_detector_nodes();

    std::set<GraphFillRegion *> regions = s.find_all_regions();

    for (auto *r : regions) {
        const char *color;
        if (r->blossom_parent_top->radius.is_growing()) {
            color = "#FF000080";
        } else if (r->blossom_parent_top->radius.is_frozen()) {
            if (r->blossom_parent_top->match.edge.loc_to == nullptr) {
                color = "#00200080";
            } else {
                color = "#00FF0080";
            }
        } else {
            color = "#0000FF80";
        }
        s.fill_region(r, color, "none");
    }

    s.draw_detection_events();
    s.draw_pending_collisions();

    s.draw_alt_tree_edges();
    s.draw_blossom_edges();
    s.draw_match_edges();

    if (focused_event.event_type == REGION_HIT_REGION) {
        auto &e = focused_event.region_hit_region_event_data.edge;
        size_t k1 = e.loc_from - &s.ns[0];
        size_t k2 = e.loc_to - &s.ns[0];
        out << " <line x1=\"" << coords[k1].first * s.scale
            << "\" x2=\"" << coords[k2].first * s.scale
            << "\" y1=\"" << coords[k1].second * s.scale
            << "\" y2=\"" << coords[k2].second * s.scale
            << "\" stroke-width=\"" << 9
            << "\" stroke=\"#000\"/>\n";
        out << " <line x1=\"" << coords[k1].first * s.scale
            << "\" x2=\"" << coords[k2].first * s.scale
            << "\" y1=\"" << coords[k1].second * s.scale
            << "\" y2=\"" << coords[k2].second * s.scale
            << "\" stroke-width=\"" << 6
            << "\" stroke=\"#F80\"/>\n";
    } else if (focused_event.event_type == REGION_HIT_BOUNDARY) {
        s.fill_region(focused_event.region_hit_boundary_event_data.region, "#80800080", "#000");
    } else if (focused_event.event_type == BLOSSOM_SHATTER) {
        s.fill_region(focused_event.blossom_shatter_event_data.blossom_region, "#FF800080", "#FF8000");
    }

    out << "</svg>\n";
}
