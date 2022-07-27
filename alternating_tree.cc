#include "alternating_tree.h"

#include <utility>


pm::AltTreeEdge::AltTreeEdge() : alt_tree_node(nullptr), edge(nullptr, nullptr, 0) {}

pm::AltTreeEdge::AltTreeEdge(AltTreeNode *alt_tree_node, const CompressedEdge &edge)
    : alt_tree_node(alt_tree_node), edge(edge) {}

void pm::AltTreeNode::add_child(const pm::AltTreeEdge& child) {
    children.push_back(child);
    child.alt_tree_node->parent = {this, child.edge.reversed()};
}

pm::AltTreeNode *pm::AltTreeNode::make_child(pm::GraphFillRegion *inner_region, pm::GraphFillRegion *outer_region,
                                             const pm::CompressedEdge &inner_to_outer_edge,
                                             const pm::CompressedEdge &child_edge) {
    auto child = new AltTreeNode(inner_region, outer_region, inner_to_outer_edge);
    auto child_alt_treeedge = AltTreeEdge(child, child_edge);
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
    outer_region->alt_tree_node = nullptr;
    inner_region->alt_tree_node = nullptr;
}

