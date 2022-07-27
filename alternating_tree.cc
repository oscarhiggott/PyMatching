#include "alternating_tree.h"

#include <utility>
#include <algorithm>


pm::AltTreeEdge::AltTreeEdge() : alt_tree_node(nullptr), edge(nullptr, nullptr, 0) {}

pm::AltTreeEdge::AltTreeEdge(AltTreeNode *alt_tree_node, const CompressedEdge &edge)
    : alt_tree_node(alt_tree_node), edge(edge) {}

bool pm::AltTreeEdge::operator==(const pm::AltTreeEdge &rhs) const {
    return alt_tree_node == rhs.alt_tree_node &&
           edge == rhs.edge;
}

bool pm::AltTreeEdge::operator!=(const pm::AltTreeEdge &rhs) const {
    return !(rhs == *this);
}

void pm::AltTreeNode::add_child(const pm::AltTreeEdge& child) {
    children.push_back(child);
    child.alt_tree_node->parent = {this, child.edge.reversed()};
}

pm::AltTreeNode *pm::AltTreeNode::make_child(pm::GraphFillRegion *child_inner_region, pm::GraphFillRegion *child_outer_region,
                                             const pm::CompressedEdge &child_inner_to_outer_edge,
                                             const pm::CompressedEdge &child_compressed_edge) {
    auto child = new AltTreeNode(child_inner_region, child_outer_region,
                                 child_inner_to_outer_edge);
    auto child_alt_treeedge = AltTreeEdge(child, child_compressed_edge);
    add_child(child_alt_treeedge);
    return child;
}

pm::AltTreeNode::AltTreeNode()
    : inner_region(nullptr), outer_region(nullptr) {}

pm::AltTreeNode::AltTreeNode(pm::GraphFillRegion *inner_region, pm::GraphFillRegion *outer_region,
                             const pm::CompressedEdge &inner_to_outer_edge,
                             const pm::AltTreeEdge &parent, std::vector<AltTreeEdge> children
                             )
      : inner_region(inner_region), outer_region(outer_region), inner_to_outer_edge(inner_to_outer_edge),
      parent(parent), children(std::move(children)){
    inner_region->alt_tree_node = this;
    outer_region->alt_tree_node = this;
}

pm::AltTreeNode::AltTreeNode(pm::GraphFillRegion *inner_region, pm::GraphFillRegion *outer_region,
                             const pm::CompressedEdge &inner_to_outer_edge)
                             : inner_region(inner_region), outer_region(outer_region),
                             inner_to_outer_edge(inner_to_outer_edge) {
    inner_region->alt_tree_node = this;
    outer_region->alt_tree_node = this;
}

pm::AltTreeNode::AltTreeNode(pm::GraphFillRegion *outer_region)
        : inner_region(nullptr), outer_region(outer_region) {}

pm::AltTreeNode *pm::AltTreeNode::find_root() {
    auto current = this;
    while (current->parent.alt_tree_node) {
        current = current->parent.alt_tree_node;
    }
    return current;
}

pm::AltTreeNode::~AltTreeNode() {
    if (outer_region)
        outer_region->alt_tree_node = nullptr;
    if (inner_region)
        inner_region->alt_tree_node = nullptr;
}


void pm::AltTreeNode::become_root() {
    // Performs a tree rotation turning this node into the root of the tree/
    auto old_parent = parent.alt_tree_node;
    old_parent->become_root();
    parent.alt_tree_node->inner_region = inner_region;
    parent.alt_tree_node->inner_to_outer_edge = parent.edge;
    inner_region = nullptr;
//    parent.alt_tree_node->children.
}

bool pm::AltTreeNode::operator==(const pm::AltTreeNode &rhs) const {
    return inner_region == rhs.inner_region &&
           outer_region == rhs.outer_region &&
           inner_to_outer_edge == rhs.inner_to_outer_edge &&
           parent == rhs.parent &&
           children == rhs.children;
}

bool pm::AltTreeNode::operator!=(const pm::AltTreeNode &rhs) const {
    return !(rhs == *this);
}

