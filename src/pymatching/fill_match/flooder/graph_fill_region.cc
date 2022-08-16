#include "pymatching/fill_match/flooder/graph_fill_region.h"

#include "pymatching/fill_match/flooder/graph.h"
#include "pymatching/fill_match/flooder_matcher_interop/varying.h"

using namespace pm;

GraphFillRegion::GraphFillRegion()
    : blossom_parent(nullptr), alt_tree_node(nullptr), radius((0 << 2) + 1), shrink_event_tracker() {
}

GraphFillRegion *GraphFillRegion::top_region() const {
    auto current = const_cast<GraphFillRegion *>(this);
    while (current->blossom_parent) {
        current = current->blossom_parent;
    }
    return current;
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
