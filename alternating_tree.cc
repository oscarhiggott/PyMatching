#include "alternating_tree.h"

#include <utility>
#include <algorithm>


pm::AltTreeEdge::AltTreeEdge() : alt_tree_node(nullptr), edge(nullptr, nullptr, 0) {}

pm::AltTreeEdge::AltTreeEdge(AltTreeNode *alt_tree_node, const CompressedEdge &edge)
    : alt_tree_node(alt_tree_node), edge(edge) {}

bool pm::AltTreeEdge::operator==(const pm::AltTreeEdge &rhs) const {
    return edge == rhs.edge;
}

bool pm::AltTreeEdge::operator!=(const pm::AltTreeEdge &rhs) const {
    return !(rhs == *this);
}

void pm::AltTreeNode::add_child(const pm::AltTreeEdge& child) {
    children.push_back(child);
    child.alt_tree_node->parent = {this, child.edge.reversed()};
}

pm::AltTreeNode *pm::AltTreeNode::make_child(pm::GraphFillRegion *child_inner_region,
                                             pm::GraphFillRegion *child_outer_region,
                                             const pm::CompressedEdge &child_inner_to_outer_edge,
                                             const pm::CompressedEdge &child_compressed_edge) {
    auto child = new AltTreeNode(child_inner_region, child_outer_region,
                                 child_inner_to_outer_edge);
    auto child_alt_tree_edge = AltTreeEdge(child, child_compressed_edge);
    add_child(child_alt_tree_edge);
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

const pm::AltTreeNode *pm::AltTreeNode::find_root() const {
    const pm::AltTreeNode* current = this;
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
    if (!parent.alt_tree_node)
        return;
    auto old_parent = parent.alt_tree_node;
    old_parent->become_root();
    parent.alt_tree_node->inner_region = inner_region;
    parent.alt_tree_node->inner_to_outer_edge = parent.edge;
    inner_region = nullptr;
    pm::unstable_erase(parent.alt_tree_node->children,
                       [&](pm::AltTreeEdge x){return x.alt_tree_node == this;}
                       );
    parent = pm::AltTreeEdge();
    add_child(pm::AltTreeEdge(old_parent, inner_to_outer_edge));
}

bool pm::AltTreeNode::operator==(const pm::AltTreeNode &rhs) const {
    return find_root()->tree_equal(*rhs.find_root());
}

bool pm::AltTreeNode::operator!=(const pm::AltTreeNode &rhs) const {
    return !(rhs == *this);
}

bool pm::AltTreeNode::tree_equal(const pm::AltTreeNode& other) const {
    if (inner_region != other.inner_region || outer_region != other.outer_region ||
        inner_to_outer_edge != other.inner_to_outer_edge ||
        children.size() != other.children.size()){
        return false;
    }
    if (children.empty())
        return true;
    for (size_t i = 0; i < children.size(); i++) {
        if (children[i].edge != other.children[i].edge)
            return false;
        if (!children[i].alt_tree_node->tree_equal(*other.children[i].alt_tree_node))
            return false;
    }
    return true;
}

std::vector<pm::AltTreeNode *> pm::AltTreeNode::all_nodes_in_tree() {
    std::vector<pm::AltTreeNode*> all_nodes;
    std::vector<pm::AltTreeNode*> to_visit;
    to_visit.push_back(this);
    while (!to_visit.empty()) {
        pm::AltTreeNode* current = to_visit.back();
        to_visit.pop_back();
        all_nodes.push_back(current);
        for (auto c: current->children)
            to_visit.push_back(c.alt_tree_node);
    }
    return all_nodes;
}
