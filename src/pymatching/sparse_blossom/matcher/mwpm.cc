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

#include "pymatching/sparse_blossom/matcher/mwpm.h"

#include <set>

#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"
#include "pymatching/sparse_blossom/matcher/alternating_tree.h"

using namespace pm;

Mwpm::Mwpm(GraphFlooder flooder) : flooder(std::move(flooder)) {
}

Mwpm::Mwpm(GraphFlooder flooder, SearchFlooder search_flooder)
    : flooder(std::move(flooder)), search_flooder(std::move(search_flooder)) {
}

Mwpm &Mwpm::operator=(Mwpm &&other) noexcept {
    if (this != &other) {
        this->~Mwpm();
        new (this) Mwpm(std::move(other));
    }
    return *this;
}

Mwpm::Mwpm(Mwpm &&other) noexcept
    : flooder(std::move(other.flooder)),
      node_arena(std::move(other.node_arena)),
      search_flooder(std::move(other.search_flooder)) {
}

void Mwpm::shatter_descendants_into_matches_and_freeze(AltTreeNode &alt_tree_node) {
    for (auto &child_edge : alt_tree_node.children) {
        shatter_descendants_into_matches_and_freeze(*child_edge.alt_tree_node);
    }
    if (alt_tree_node.inner_region) {
        alt_tree_node.parent = AltTreeEdge();
        alt_tree_node.inner_region->add_match(alt_tree_node.outer_region, alt_tree_node.inner_to_outer_edge);
        flooder.set_region_frozen(*alt_tree_node.inner_region);
        flooder.set_region_frozen(*alt_tree_node.outer_region);
        alt_tree_node.inner_region->alt_tree_node = nullptr;
        alt_tree_node.outer_region->alt_tree_node = nullptr;
    }
    if (alt_tree_node.outer_region) {
        alt_tree_node.outer_region->alt_tree_node = nullptr;
    }
    node_arena.del(&alt_tree_node);
}

void Mwpm::handle_tree_hitting_boundary(const RegionHitBoundaryEventData &event) {
    auto node = event.region->alt_tree_node;
    node->become_root();
    // Match descendents, deleting AltTreeNodes and freezing GraphFillRegions
    shatter_descendants_into_matches_and_freeze(*node);

    // Now match the event region to the boundary and freeze
    event.region->match = Match{nullptr, event.edge};
    flooder.set_region_frozen(*event.region);
}

void Mwpm::handle_tree_hitting_boundary_match(
    GraphFillRegion *unmatched_region,
    GraphFillRegion *matched_region,
    const CompressedEdge &unmatched_to_matched_edge) {
    auto &alt_tree_node = unmatched_region->alt_tree_node;
    unmatched_region->add_match(matched_region, unmatched_to_matched_edge);
    flooder.set_region_frozen(*unmatched_region);
    alt_tree_node->become_root();
    shatter_descendants_into_matches_and_freeze(*alt_tree_node);
}

void Mwpm::handle_tree_hitting_other_tree(const RegionHitRegionEventData &event) {
    auto alt_node_1 = event.region1->alt_tree_node;
    auto alt_node_2 = event.region2->alt_tree_node;
    // Tree rotation
    event.region1->alt_tree_node->become_root();
    event.region2->alt_tree_node->become_root();
    // Match and freeze descendants
    shatter_descendants_into_matches_and_freeze(*alt_node_1);
    shatter_descendants_into_matches_and_freeze(*alt_node_2);
    // Match colliding nodes
    event.region1->add_match(event.region2, event.edge);
    // Freeze colliding regions
    flooder.set_region_frozen(*event.region1);
    flooder.set_region_frozen(*event.region2);
}

AltTreeNode *Mwpm::make_child(
    AltTreeNode &parent,
    GraphFillRegion *child_inner_region,
    GraphFillRegion *child_outer_region,
    const CompressedEdge &child_inner_to_outer_edge,
    const CompressedEdge &child_compressed_edge) {
    auto child = node_arena.alloc_unconstructed();
    new (child) AltTreeNode(child_inner_region, child_outer_region, child_inner_to_outer_edge);
    auto child_alt_tree_edge = AltTreeEdge(child, child_compressed_edge);
    parent.add_child(child_alt_tree_edge);
    return child;
}

void Mwpm::handle_tree_hitting_match(
    GraphFillRegion *unmatched_region,
    GraphFillRegion *matched_region,
    const CompressedEdge &unmatched_to_matched_edge) {
    auto alt_tree_node = unmatched_region->alt_tree_node;
    make_child(
        *alt_tree_node,
        matched_region,
        matched_region->match.region,
        matched_region->match.edge,
        unmatched_to_matched_edge);
    auto other_match = matched_region->match.region;
    other_match->match.clear();
    matched_region->match.clear();
    flooder.set_region_shrinking(*matched_region);
    flooder.set_region_growing(*other_match);
}

void Mwpm::handle_tree_hitting_self(const RegionHitRegionEventData &event, AltTreeNode *common_ancestor) {
    auto alt_node_1 = event.region1->alt_tree_node;
    auto alt_node_2 = event.region2->alt_tree_node;
    auto prune_result_1 = alt_node_1->prune_upward_path_stopping_before(node_arena, common_ancestor, true);
    auto prune_result_2 = alt_node_2->prune_upward_path_stopping_before(node_arena, common_ancestor, false);

    // Construct blossom region cycle
    auto blossom_cycle = std::move(prune_result_2.pruned_path_region_edges);
    auto p1s = prune_result_1.pruned_path_region_edges.size();
    blossom_cycle.reserve(blossom_cycle.size() + p1s + 1);
    for (size_t i = 0; i < p1s; i++)
        blossom_cycle.push_back(prune_result_1.pruned_path_region_edges[p1s - i - 1]);
    blossom_cycle.push_back(RegionEdge{event.region1, event.edge});
    common_ancestor->outer_region->alt_tree_node = nullptr;
    auto blossom_region = flooder.create_blossom(blossom_cycle);

    common_ancestor->outer_region = blossom_region;
    blossom_region->alt_tree_node = common_ancestor;
    common_ancestor->children.reserve(
        common_ancestor->children.size() + prune_result_1.orphan_edges.size() + prune_result_2.orphan_edges.size());
    for (auto &c : prune_result_1.orphan_edges) {
        common_ancestor->add_child(c);
    }
    for (auto &c : prune_result_2.orphan_edges) {
        common_ancestor->add_child(c);
    }
}

void Mwpm::handle_blossom_shattering(const BlossomShatterEventData &event) {
    for (auto &child : event.blossom_region->blossom_children) {
        child.region->clear_blossom_parent();
    }

    // First find indices of in_parent_region and in_child_region
    // in_parent_region is the blossom cycle region connected to the parent of the blossom inner node.
    // in_child_region is the blossom cycle region connected to the child of the inner node
    auto blossom_cycle = std::move(event.blossom_region->blossom_children);
    auto blossom_alt_node = event.blossom_region->alt_tree_node;
    size_t bsize = blossom_cycle.size();
    size_t parent_idx = 0;
    size_t child_idx = 0;
    for (size_t i = 0; i < bsize; i++) {
        if (blossom_cycle[i].region == event.in_parent_region) {
            parent_idx = i;
        }
        if (blossom_cycle[i].region == event.in_child_region) {
            child_idx = i;
        }
    }

    // Length of path starting on in_parent and stopping before in_child
    size_t gap = ((child_idx + bsize - parent_idx) % bsize);
    AltTreeNode *current_alt_node;
    size_t evens_start, evens_end;

    current_alt_node = event.blossom_region->alt_tree_node->parent.alt_tree_node;
    unstable_erase(current_alt_node->children, [blossom_alt_node](AltTreeEdge x) {
        return x.alt_tree_node == blossom_alt_node;
    });
    auto child_edge = blossom_alt_node->parent.edge.reversed();

    if (gap % 2 == 0) {
        // The path starting after in_child and stopping before in_parent is even length. Regions will
        // be matched along this path
        evens_start = child_idx + 1;
        evens_end = child_idx + bsize - gap;

        // Now insert odd-length path starting on in_parent and ending on in_child into alternating tree
        for (size_t i = parent_idx; i < parent_idx + gap; i += 2) {
            current_alt_node = make_child(
                *current_alt_node,
                blossom_cycle[i % bsize].region,
                blossom_cycle[(i + 1) % bsize].region,
                blossom_cycle[i % bsize].edge,
                child_edge);
            child_edge = blossom_cycle[(i + 1) % bsize].edge;
            flooder.set_region_shrinking(*current_alt_node->inner_region);
            flooder.set_region_growing(*current_alt_node->outer_region);
        }
    } else {
        // The path starting after in_parent and stopping before in_child is even length. Regions will
        // be matched along this path.
        evens_start = parent_idx + 1;
        evens_end = parent_idx + gap;

        // Now insert odd-length path into alternating tree
        for (size_t i = 0; i < bsize - gap; i += 2) {
            size_t k1 = (parent_idx + bsize - i) % bsize;
            size_t k2 = (parent_idx + bsize - i - 1) % bsize;
            size_t k3 = (parent_idx + bsize - i - 2) % bsize;
            current_alt_node = make_child(
                *current_alt_node,
                blossom_cycle[k1].region,
                blossom_cycle[k2].region,
                blossom_cycle[k2].edge.reversed(),
                child_edge);
            child_edge = blossom_cycle[k3].edge.reversed();
            flooder.set_region_shrinking(*current_alt_node->inner_region);
            flooder.set_region_growing(*current_alt_node->outer_region);
        }
    }

    for (size_t j = evens_start; j < evens_end; j += 2) {
        size_t k1 = j % bsize;
        size_t k2 = (j + 1) % bsize;
        blossom_cycle[k1].region->add_match(blossom_cycle[k2].region, blossom_cycle[k1].edge);

        // The blossom regions were previously shrinking. Now they are stopped. This can create new
        // events on the nodes, and so the nodes must be reprocessed.
        blossom_cycle[k1].region->do_op_for_each_node_in_total_area([this](DetectorNode *n) {
            flooder.reschedule_events_at_detector_node(*n);
        });
        blossom_cycle[k2].region->do_op_for_each_node_in_total_area([this](DetectorNode *n) {
            flooder.reschedule_events_at_detector_node(*n);
        });
    }

    blossom_alt_node->inner_region = blossom_cycle[child_idx].region;
    flooder.set_region_shrinking(*blossom_alt_node->inner_region);
    blossom_cycle[child_idx].region->alt_tree_node = blossom_alt_node;
    current_alt_node->add_child(AltTreeEdge(blossom_alt_node, child_edge));

    flooder.region_arena.del(event.blossom_region);
}

void Mwpm::handle_region_hit_region(const MwpmEvent event) {
    const auto &d = event.region_hit_region_event_data;
    auto alt_node_1 = d.region1->alt_tree_node;
    auto alt_node_2 = d.region2->alt_tree_node;
    if (alt_node_1 && alt_node_2) {
        auto common_ancestor = alt_node_1->most_recent_common_ancestor(*alt_node_2);
        if (!common_ancestor) {
            handle_tree_hitting_other_tree(d);
        } else {
            handle_tree_hitting_self(d, common_ancestor);
        }
    } else if (alt_node_1) {
        // Region 2 is not in the tree, so must be matched to the boundary or another region
        if (d.region2->match.region) {
            handle_tree_hitting_match(d.region1, d.region2, d.edge);
        } else {
            handle_tree_hitting_boundary_match(d.region1, d.region2, d.edge);
        }
    } else {
        // Region 1 is not in the tree, so must be matched to the boundary or another region
        if (d.region1->match.region) {
            handle_tree_hitting_match(d.region2, d.region1, d.edge.reversed());
        } else {
            handle_tree_hitting_boundary_match(d.region2, d.region1, d.edge.reversed());
        }
    }
}
void Mwpm::process_event(const MwpmEvent &event) {
    switch (event.event_type) {
        case REGION_HIT_REGION:
            handle_region_hit_region(event);
            break;
        case REGION_HIT_BOUNDARY:
            handle_tree_hitting_boundary(event.region_hit_boundary_event_data);
            break;
        case BLOSSOM_SHATTER:
            handle_blossom_shattering(event.blossom_shatter_event_data);
            break;
        case NO_EVENT:
            // Do nothing.
            break;
        default:
            throw std::invalid_argument("Unrecognized event type");
    }
}

GraphFillRegion *Mwpm::pair_and_shatter_subblossoms_and_extract_matches(GraphFillRegion *region, MatchingResult &res) {
    for (auto &r : region->blossom_children) {
        r.region->clear_blossom_parent_ignoring_wrapped_radius();
    }
    auto subblossom = region->match.edge.loc_from->region_that_arrived_top;
    subblossom->match = region->match;
    if (subblossom->match.region)
        subblossom->match.region->match.region = subblossom;
    res.weight += region->radius.y_intercept();
    auto iter = std::find_if(
        region->blossom_children.begin(), region->blossom_children.end(), [&subblossom](const RegionEdge &e) {
            return e.region == subblossom;
        });
    size_t index = std::distance(region->blossom_children.begin(), iter);
    size_t num_children = region->blossom_children.size();
    for (size_t i = 0; i < num_children - 1; i += 2) {
        auto &re1 = region->blossom_children[(index + i + 1) % num_children];
        auto &re2 = region->blossom_children[(index + i + 2) % num_children];
        re1.region->add_match(re2.region, re1.edge);
        res += shatter_blossom_and_extract_matches(re1.region);
    }
    flooder.region_arena.del(region);
    return subblossom;
}

MatchingResult Mwpm::shatter_blossom_and_extract_matches(GraphFillRegion *region) {
    region->cleanup_shell_area();

    // First handle base cases (no subblossoms)
    if (region->match.region) {
        region->match.region->cleanup_shell_area();
        if (region->blossom_children.empty() && region->match.region->blossom_children.empty()) {
            // Neither region nor matched region have blossom children
            // No shattering required, so just return MatchingResult from this match.
            MatchingResult res = {
                region->match.edge.obs_mask, region->radius.y_intercept() + region->match.region->radius.y_intercept()};
            flooder.region_arena.del(region->match.region);
            flooder.region_arena.del(region);
            return res;
        }
    } else if (region->blossom_children.empty()) {
        // Region with no blossom children matched to boundary
        // No shattering required, so just return MatchingResult from this match.
        MatchingResult res = {region->match.edge.obs_mask, region->radius.y_intercept()};
        flooder.region_arena.del(region);
        return res;
    }

    // Pair up and shatter subblossoms into matches
    MatchingResult res{0, 0};
    if (!region->blossom_children.empty())
        region = pair_and_shatter_subblossoms_and_extract_matches(region, res);
    if (region->match.region && !region->match.region->blossom_children.empty())
        pair_and_shatter_subblossoms_and_extract_matches(region->match.region, res);
    res += shatter_blossom_and_extract_matches(region);
    return res;
}

GraphFillRegion *Mwpm::pair_and_shatter_subblossoms_and_extract_match_edges(
    GraphFillRegion *region, std::vector<CompressedEdge> &match_edges) {
    for (auto &r : region->blossom_children) {
        r.region->clear_blossom_parent_ignoring_wrapped_radius();
    }
    auto subblossom = region->match.edge.loc_from->region_that_arrived_top;
    subblossom->match = region->match;
    if (subblossom->match.region)
        subblossom->match.region->match.region = subblossom;
    auto iter = std::find_if(
        region->blossom_children.begin(), region->blossom_children.end(), [&subblossom](const RegionEdge &e) {
            return e.region == subblossom;
        });
    size_t index = std::distance(region->blossom_children.begin(), iter);
    size_t num_children = region->blossom_children.size();
    for (size_t i = 0; i < num_children - 1; i += 2) {
        auto &re1 = region->blossom_children[(index + i + 1) % num_children];
        auto &re2 = region->blossom_children[(index + i + 2) % num_children];
        re1.region->add_match(re2.region, re1.edge);
        shatter_blossom_and_extract_match_edges(re1.region, match_edges);
    }
    flooder.region_arena.del(region);
    return subblossom;
}

void Mwpm::shatter_blossom_and_extract_match_edges(GraphFillRegion *region, std::vector<CompressedEdge> &match_edges) {
    region->cleanup_shell_area();

    // First handle base cases (no subblossoms)
    if (region->match.region) {
        region->match.region->cleanup_shell_area();
        if (region->blossom_children.empty() && region->match.region->blossom_children.empty()) {
            // Neither region nor matched region have blossom children
            // No shattering required, so just return MatchingResult from this match.
            match_edges.push_back(region->match.edge);
            flooder.region_arena.del(region->match.region);
            flooder.region_arena.del(region);
            return;
        }
    } else if (region->blossom_children.empty()) {
        // Region with no blossom children matched to boundary
        // No shattering required, so just return MatchingResult from this match.
        match_edges.push_back(region->match.edge);
        flooder.region_arena.del(region);
        return;
    }

    // Pair up and shatter subblossoms into matches
    if (!region->blossom_children.empty())
        region = pair_and_shatter_subblossoms_and_extract_match_edges(region, match_edges);
    if (region->match.region && !region->match.region->blossom_children.empty())
        pair_and_shatter_subblossoms_and_extract_match_edges(region->match.region, match_edges);
    shatter_blossom_and_extract_match_edges(region, match_edges);
    return;
}

void Mwpm::create_detection_event(DetectorNode *node) {
    auto region = flooder.region_arena.alloc_default_constructed();
    auto alt_tree_node = node_arena.alloc_unconstructed();
    new (alt_tree_node) AltTreeNode(region);
    region->alt_tree_node = alt_tree_node;
    flooder.do_region_created_at_empty_detector_node(*region, *node);
}

MatchingResult &MatchingResult::operator+=(const MatchingResult &rhs) {
    obs_mask ^= rhs.obs_mask;
    weight += rhs.weight;
    return *this;
}

MatchingResult MatchingResult::operator+(const MatchingResult &rhs) const {
    MatchingResult copy = *this;
    copy += rhs;
    return copy;
}

MatchingResult::MatchingResult() : obs_mask(0), weight(0) {
}

MatchingResult::MatchingResult(obs_int obs_mask, total_weight_int weight) : obs_mask(obs_mask), weight(weight) {
}

bool MatchingResult::operator==(const MatchingResult &rhs) const {
    return obs_mask == rhs.obs_mask && weight == rhs.weight;
}

bool MatchingResult::operator!=(const MatchingResult &rhs) const {
    return !(rhs == *this);
}

void Mwpm::verify_invariants() const {
}

void Mwpm::extract_paths_from_match_edges(
    const std::vector<CompressedEdge> &match_edges, uint8_t *obs_begin_ptr, total_weight_int &weight) {
    for (auto &edge : match_edges) {
        size_t loc_to_idx = edge.loc_to ? edge.loc_to - &flooder.graph.nodes[0] : SIZE_MAX;
        search_flooder.iter_edges_on_shortest_path_from_middle(
            edge.loc_from - &flooder.graph.nodes[0], loc_to_idx, [&](const pm::SearchGraphEdge &e) {
                auto &obs = e.detector_node->neighbor_observable_indices[e.neighbor_index];
                for (auto i : obs)
                    *(obs_begin_ptr + i) ^= 1;
                weight += e.detector_node->neighbor_weights[e.neighbor_index];
            });
    }
}

Mwpm::Mwpm() {
}

void Mwpm::reset() {
    for (auto &n : flooder.graph.nodes)
        n.reset();
    for (auto &m : search_flooder.graph.nodes)
        m.reset();
    flooder.queue.clear();
    node_arena.~Arena();
    flooder.region_arena.~Arena();
}
