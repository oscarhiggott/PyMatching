#ifndef PYMATCHING2_GRAPH_FILL_REGION_H
#define PYMATCHING2_GRAPH_FILL_REGION_H

#include <vector>

#include "pymatching/fill_match/flooder_matcher_interop/varying.h"
#include "pymatching/fill_match/flooder_matcher_interop/region_edge.h"
#include "pymatching/fill_match/flooder/match.h"
#include "pymatching/fill_match/tracker/queued_event_tracker.h"

namespace pm {

class AltTreeNode;

struct GraphFillRegion {
    /// If this region has merged with others into a blossom, this is that blossom.
    GraphFillRegion* blossom_parent;
    /// The topmost fill region that contains this region. This field must be kept up to date as
    /// the region structure changes.
    GraphFillRegion* blossom_parent_top;
    /// If this is a top-level region (not a blossom child), this is the alternating tree that
    /// it is part of. Note that it may be a degenerate alternating tree with just a single
    /// graph fill region (this one).
    pm::AltTreeNode* alt_tree_node;
    /// How much this region has grown since its creation (and whether it is currently growing).
    /// For graph fill regions starting from detection events, this is just the actual radius
    /// of the region. For graph fill regions starting from blossom-creation events, this is how
    /// much the blossom has grown since it was created (it's the *extra* radius, not the total
    /// radius starting from the detection events).
    pm::Varying32 radius;
    /// Event tracker for shrink events. Handles ensuring at least, and ideally exactly, one event
    /// to look at the region is in the event queue.
    QueuedEventTracker shrink_event_tracker;
    /// If the region is matched, as opposed to growing/shrinking, this says what it is matched to.
    pm::Match match;

    /// If this region is a blossom, these are its child regions along with the cyclic paths
    /// between the children.
    std::vector<pm::RegionEdge> blossom_children;
    /// The set of nodes directly owned by this region. Note that this vector does not include
    /// the nodes indirectly owned by this region that are owned by the blossom children of this
    /// region (or their children or etc).
    std::vector<pm::DetectorNode*> shell_area;

    void cleanup_shell_area();

    GraphFillRegion();
    GraphFillRegion(GraphFillRegion &&);
    GraphFillRegion(const GraphFillRegion &) = delete;
    bool tree_equal(const pm::GraphFillRegion& other) const;

    void add_match(pm::GraphFillRegion* match, const pm::CompressedEdge& edge);

    template <typename Callable>
    void do_op_for_each_node_in_total_area(const Callable& func);
    template <typename Callable>
    void do_op_for_each_descendant_and_self(const Callable& func);
    void set_blossom_parent(GraphFillRegion *new_blossom_parent);

    bool operator==(const GraphFillRegion& rhs) const;

    bool operator!=(const GraphFillRegion& rhs) const;
};

template <typename Callable>
inline void pm::GraphFillRegion::do_op_for_each_node_in_total_area(const Callable& func) {
    for (size_t i = 0; i < shell_area.size(); i++) {
        func(shell_area[shell_area.size() - i - 1]);
    }
    for (auto& child : blossom_children) {
        child.region->do_op_for_each_node_in_total_area(func);
    }
}

template <typename Callable>
inline void pm::GraphFillRegion::do_op_for_each_descendant_and_self(const Callable& func) {
    func(this);
    for (RegionEdge &child : blossom_children) {
        child.region->do_op_for_each_descendant_and_self(func);
    }
}

}  // namespace pm

#endif  // PYMATCHING2_GRAPH_FILL_REGION_H
