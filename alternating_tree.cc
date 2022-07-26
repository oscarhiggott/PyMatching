#include "alternating_tree.h"


pm::AltTreeEdge::AltTreeEdge() : alt_tree_node(nullptr), edge(nullptr, nullptr, 0) {}

pm::AltTreeEdge::AltTreeEdge(AltTreeNode *alt_tree_node, const CompressedEdge &edge)
    : alt_tree_node(alt_tree_node), edge(edge) {}

void pm::AltTreeNode::add_child(const pm::AltTreeEdge& child) {
    children.push_back(child);
    child.alt_tree_node->parent = {this, child.edge.reversed()};
}

pm::AltTreeNode *pm::AltTreeNode::make_child(pm::GraphFillRegion *inner_region, pm::GraphFillRegion *outer_region,
                                             const pm::CompressedEdge &inner_outer_edge,
                                             const pm::CompressedEdge &child_edge) {
    return nullptr;
}

pm::AltTreeNode::AltTreeNode()
    : inner_region(nullptr), outer_region(nullptr) {

}

pm::AltTreeNode::AltTreeNode(pm::GraphFillRegion *inner_region, pm::GraphFillRegion *outer_region,
                             const pm::AltTreeEdge &parent, const std::vector<AltTreeEdge> &children,
                             const pm::CompressedEdge &inner_to_outer_edge)
      : inner_region(inner_region), outer_region(outer_region), parent(parent), children(children),
      inner_to_outer_edge(inner_to_outer_edge) {}
