#include "alternating_tree.h"


pm::AltTreeEdge::AltTreeEdge() : alt_tree_node(nullptr), edge(nullptr, nullptr, 0) {}

pm::AltTreeEdge::AltTreeEdge(AltTreeNode *alt_tree_node, const CompressedEdge &edge)
    : alt_tree_node(alt_tree_node), edge(edge) {}

void pm::AltTreeNode::add_child(const pm::AltTreeEdge& child) {
    children.push_back(child);
    child.alt_tree_node->parent = {this, child.edge.reversed()};
}
