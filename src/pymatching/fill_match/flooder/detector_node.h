#ifndef PYMATCHING_FILL_MATCH_DETECTOR_NODE_H
#define PYMATCHING_FILL_MATCH_DETECTOR_NODE_H

#include <vector>

#include "pymatching/fill_match/flooder_matcher_interop/varying.h"
#include "pymatching/fill_match/tracker/queued_event_tracker.h"
#include "pymatching/fill_match/flooder/graph_fill_region.h"

namespace pm {

/// A detector node is a location where a detection event might occur.
///
/// It corresponds to a potential symptom that could be seen, and can
/// be used to infer which errors have occurred.
///
/// There is a 1:1 correspondence between detector nodes in a matching
/// graph and the DETECTOR annotations in a Stim circuit.
class DetectorNode {
   public:
    DetectorNode()
        : observables_crossed_from_source(0),
          reached_from_source(nullptr),
          radius_of_arrival(0),
          region_that_arrived(nullptr),
          region_that_arrived_top(nullptr),
          wrapped_radius_cached(0) {
    }

    /// == Ephemeral fields used to track algorithmic state during matching. ==
    /// The region that reached and owns this node.
    GraphFillRegion* region_that_arrived;
    /// The topmost region containing this node. Must be kept up to date as
    /// the region structure changes.
    GraphFillRegion* region_that_arrived_top;
    /// Stores the latest value of `compute_wrapped_radius` so it doesn't need to be recomputed
    /// as much. Updated whenever region_that_arrived_top changes.
    int32_t wrapped_radius_cached;
    /// Of the detection events within the owning region, which one is this node linked to.
    DetectorNode* reached_from_source;
    /// Which observables are crossed, travelling from this node to the source
    /// detection event that reached it. Must be 0 if reached_from_source == nullptr.
    obs_int observables_crossed_from_source;
    /// This was the radius of `region_that_arrived` when it arrived at this node. Otherwise 0.
    cumulative_time_int radius_of_arrival;
    /// Manages the next "look at me!" event for the node.
    QueuedEventTracker node_event_tracker;

    /// == Permanent fields used to define the structure of the graph. ==
    std::vector<DetectorNode*> neighbors;       /// The node's neighbors.
    std::vector<weight_int> neighbor_weights;   /// Distance crossed by the edge to each neighbor.
    std::vector<obs_int> neighbor_observables;  /// Observables crossed by the edge to each neighbor.

    /// After it reached this node, how much further did the owning search region grow? Also is it currently growing?
    inline Varying32 local_radius() const {
        if (region_that_arrived_top == nullptr) {
            return Varying32{0};
        }
        return region_that_arrived_top->radius + wrapped_radius_cached;
    }

    /// Check if this node is part the same top-level region as another.
    /// Note that they may have different lower level owners that still merge into the same top level owned.
    inline bool has_same_owner_as(const DetectorNode &other) const {
        return region_that_arrived_top == other.region_that_arrived_top;
    }

    /// Zero all ephemeral fields.
    /// Doesn't free anything or propagate a signal to other objects. Just zeros the fields.
    void reset();

    size_t index_of_neighbor(DetectorNode* neighbor) const;

    int32_t compute_wrapped_radius() const;
};

}  // namespace pm

#endif  // PYMATCHING_FILL_MATCH_DETECTOR_NODE_H
