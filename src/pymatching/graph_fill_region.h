#ifndef PYMATCHING2_GRAPH_FILL_REGION_H
#define PYMATCHING2_GRAPH_FILL_REGION_H

#include "pymatching/alternating_tree.h"
#include "pymatching/region_edge.h"

namespace pm {

class AltTreeNode;

/// A Match is a partnered graph fill region.
///
/// When one graph fill region is matched to another, it's guaranteed that (if the matching
/// process terminated) two of their detection events will be matched to each other. The
/// specific detection events to match are identified by `edge.loc_from` and `edge.loc_to`.
struct Match {
    /// The region being matched to (or nullptr if not matched).
    pm::GraphFillRegion* region;
    /// A summary of the low-level path from this region to that region, including the
    /// start/end detection events and the observables that were crossed.
    pm::CompressedEdge edge;

    Match();
    Match(pm::GraphFillRegion* region, pm::CompressedEdge edge);

    bool operator==(const Match& rhs) const;
    bool operator!=(const Match& rhs) const;
};

class GraphFillRegion {
   public:
    // TODO: Maybe one field for blossom_parent or alt_tree_node eventually by using union?

    /// If this region has merged with others into a blossom, this is that blossom.
    GraphFillRegion* blossom_parent;
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
    /// Event validation index for blossom shatter event.
    uint64_t shrink_event_vid;
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
    bool tree_equal(const pm::GraphFillRegion& other) const;

    GraphFillRegion* top_region() const;
    void add_match(pm::GraphFillRegion* match, const pm::CompressedEdge& edge);

    template <typename Callable>
    void do_op_for_each_node_in_total_area(const Callable& func);

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

}  // namespace pm

#endif  // PYMATCHING2_GRAPH_FILL_REGION_H
