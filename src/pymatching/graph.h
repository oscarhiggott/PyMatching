#ifndef PYMATCHING2_GRAPH_H
#define PYMATCHING2_GRAPH_H

#include <vector>

#include "pymatching/bit_bucket_queue.h"
#include "pymatching/fixed_length_vector.h"
#include "pymatching/varying.h"

namespace pm {

class GraphFillRegion;
struct TentativeEvent;

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
          distance_from_source(0),
          region_that_arrived(nullptr) {
    }

    /// == Ephemeral fields used to track algorithmic state during matching. ==
    /// The region that reached and owns this node.
    GraphFillRegion* region_that_arrived;
    /// Of the detection events within the owning region, which one is this node linked to.
    DetectorNode* reached_from_source;
    /// Which observables are crossed, travelling from this node to the source
    /// detection event that reached it. Must be 0 if reached_from_source == nullptr.
    obs_int observables_crossed_from_source;
    /// How far is it from this node to the source detection event that reached it.
    /// Must be 0 if reached_from_source == nullptr.
    cumulative_time_int distance_from_source;
    /// Manages the next "look at me!" event for the node.
    QueuedEventTracker node_event_tracker;

    /// == Permanent fields used to define the structure of the graph. ==
    std::vector<DetectorNode*> neighbors;       /// The node's neighbors.
    std::vector<weight_int> neighbor_weights;   /// Distance crossed by the edge to each neighbor.
    std::vector<obs_int> neighbor_observables;  /// Observables crossed by the edge to each neighbor.

    /// After it reached this node, how much further did the owning search region grow? Also is it currently growing?
    Varying32 local_radius() const;
    /// What's the top-level region that owns this node? The blossom-with-no-blossom-parent that this node is inside.
    GraphFillRegion* top_region() const;
    /// Check if this node is part the same top-level region as another.
    /// Note that they may have different lower level owners that still merge into the same top level owned.
    bool has_same_owner_as(const DetectorNode& other) const;
    /// Zero all ephemeral fields.
    /// Doesn't free anything or propagate a signal to other objects. Just zeros the fields.
    void reset();

    size_t index_of_neighbor(DetectorNode *neighbor) const;

    /// Implementation detail of local_radius.
    ///
    /// Tallies up the distance all the way from this detector node to the detection event that
    /// reached it first. Adds up the incremental radius of each blossom region nested on top of
    /// that detection event until reaching the top-level region that contains this node.
    Varying32 total_radius() const;
};

/// A collection of detector nodes. It's expected that all detector nodes in the graph
/// will only refer to other detector nodes within the same graph.
class MatchingGraph {
   public:
    std::vector<DetectorNode> nodes;
    size_t num_nodes;

    MatchingGraph();
    explicit MatchingGraph(size_t num_nodes);
    MatchingGraph(MatchingGraph&& graph) noexcept;
    void add_edge(size_t u, size_t v, weight_int weight, obs_int observables);
    void add_boundary_edge(size_t u, weight_int weight, obs_int observables);
};

}  // namespace pm

#endif  // PYMATCHING2_GRAPH_H
