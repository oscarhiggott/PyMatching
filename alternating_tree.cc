#include "alternating_tree.h"


void pm::AltTreeNode::add_child(const pm::AltTreeEdge& child) {
    children.push_back(child);
    child.alt_tree_node->parent = {this, child.edge.reversed()};
}