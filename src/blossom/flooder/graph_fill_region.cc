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

#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"

#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"

using namespace pm;

GraphFillRegion::GraphFillRegion()
    : blossom_parent(nullptr),
      blossom_parent_top(this),
      alt_tree_node(nullptr),
      radius((0 << 2) + 1),
      shrink_event_tracker() {
}
GraphFillRegion::GraphFillRegion(GraphFillRegion &&other)
    : blossom_parent(other.blossom_parent),
      blossom_parent_top(other.blossom_parent_top == &other ? this : other.blossom_parent_top),
      alt_tree_node(std::move(other.alt_tree_node)),
      radius(std::move(other.radius)),
      shrink_event_tracker(std::move(other.shrink_event_tracker)),
      match(std::move(other.match)),
      blossom_children(std::move(other.blossom_children)),
      shell_area(std::move(other.shell_area)) {
}

bool GraphFillRegion::tree_equal(const GraphFillRegion &other) const {
    if (alt_tree_node != other.alt_tree_node || radius != other.radius ||
        blossom_children.size() != other.blossom_children.size() || shell_area != other.shell_area) {
        return false;
    }
    if (blossom_children.empty())
        return true;
    for (size_t i = 0; i < blossom_children.size(); i++) {
        if (blossom_children[i].edge != other.blossom_children[i].edge)
            return false;
        if (!blossom_children[i].region->tree_equal(*other.blossom_children[i].region))
            return false;
    }
    return true;
}

bool GraphFillRegion::operator==(const GraphFillRegion &rhs) const {
    return tree_equal(rhs);
}

bool GraphFillRegion::operator!=(const GraphFillRegion &rhs) const {
    return !(rhs == *this);
}

void GraphFillRegion::add_match(GraphFillRegion *region, const CompressedEdge &edge) {
    match = Match{region, edge};
    region->match = Match{this, edge.reversed()};
}

void GraphFillRegion::cleanup_shell_area() {
    for (auto &detector_node : shell_area) {
        detector_node->reset();
    }
}

void GraphFillRegion::clear_blossom_parent() {
    blossom_parent = nullptr;
    do_op_for_each_descendant_and_self([&](GraphFillRegion *descendant) {
        descendant->blossom_parent_top = this;
        for (DetectorNode *n : descendant->shell_area) {
            n->region_that_arrived_top = this;
            n->wrapped_radius_cached = n->compute_wrapped_radius();
        }
    });
}

void GraphFillRegion::clear_blossom_parent_ignoring_wrapped_radius() {
    blossom_parent = nullptr;
    do_op_for_each_descendant_and_self([&](GraphFillRegion *descendant) {
        descendant->blossom_parent_top = this;
        for (DetectorNode *n : descendant->shell_area) {
            n->region_that_arrived_top = this;
        }
    });
}

void GraphFillRegion::wrap_into_blossom(GraphFillRegion *new_blossom_parent_and_top) {
    blossom_parent = new_blossom_parent_and_top;
    do_op_for_each_descendant_and_self([&](GraphFillRegion *descendant) {
        descendant->blossom_parent_top = new_blossom_parent_and_top;
        for (DetectorNode *n : descendant->shell_area) {
            n->region_that_arrived_top = new_blossom_parent_and_top;
            n->wrapped_radius_cached = n->compute_wrapped_radius();
        }
    });
}

bool GraphFillRegion::operator<=(const GraphFillRegion &rhs) const {
    const GraphFillRegion *r = this;
    while (r != nullptr && r != &rhs) {
        r = r->blossom_parent;
    }
    return r == &rhs;
}

bool GraphFillRegion::operator<(const GraphFillRegion &rhs) const {
    return this != &rhs && (*this <= rhs);
}

bool GraphFillRegion::operator>(const GraphFillRegion &rhs) const {
    return rhs < *this;
}

bool GraphFillRegion::operator>=(const GraphFillRegion &rhs) const {
    return rhs <= *this;
}
