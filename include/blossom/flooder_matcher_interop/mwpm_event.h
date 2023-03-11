// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PYMATCHING_FILL_MATCH_MWPM_EVENT_H
#define PYMATCHING_FILL_MATCH_MWPM_EVENT_H

#include "pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"

namespace pm {

struct DetectorNode;
struct GraphFillRegion;

struct RegionHitRegionEventData {
    /// One of the colliding regions.
    GraphFillRegion *region1;
    /// The other colliding region.
    GraphFillRegion *region2;
    /// The path between the two regions, anchored to detection events in the regions.
    CompressedEdge edge;

    inline RegionHitRegionEventData reversed() const {
        return {region2, region1, edge.reversed()};
    }
    bool operator==(const RegionHitRegionEventData &rhs) const;
    bool operator!=(const RegionHitRegionEventData &rhs) const;
    std::string str() const;
};
std::ostream &operator<<(std::ostream &out, const RegionHitRegionEventData &dat);

struct RegionHitBoundaryEventData {
    /// The growing region that hit the boundary.
    GraphFillRegion *region;
    /// The path from the region to the boundary, anchored to a detection event in the region.
    CompressedEdge edge;

    bool operator==(const RegionHitBoundaryEventData &rhs) const;
    bool operator!=(const RegionHitBoundaryEventData &rhs) const;
    std::string str() const;
};
std::ostream &operator<<(std::ostream &out, const RegionHitBoundaryEventData &ev);

struct BlossomShatterEventData {
    /// The shrinking blossom region that has become empty and needs to be destroyed.
    GraphFillRegion *blossom_region;
    /// The blossom child that is linked to the alternating tree parent of the blossom.
    GraphFillRegion *in_parent_region;
    /// The blossom child that is linked to the alternating tree child of the blossom.
    GraphFillRegion *in_child_region;

    bool operator==(const BlossomShatterEventData &rhs) const;
    bool operator!=(const BlossomShatterEventData &rhs) const;
    std::string str() const;
};
std::ostream &operator<<(std::ostream &out, const BlossomShatterEventData &ev);

enum MwpmEventType : uint8_t { NO_EVENT, REGION_HIT_REGION, REGION_HIT_BOUNDARY, BLOSSOM_SHATTER };

/// A MwpmEvent is an interaction that the min-weight-perform-matching algorithm must react to.
///
/// An example of a tentative event that the flooder produces while running that DOES NOT produce
/// a MwpmEvent is a region growing enough to cover a new empty detector node.
///
/// An example of a tentative event tha the flooder produce while running that does produce a
/// MwpmEvent is a region growing into a detector node owned by another region, as this indicates
/// the regions have collided and a match (or blossom or etc) needs to be formed.
struct MwpmEvent {
    /// Event data. The `event_type` field determines which is populated.
    union {
        RegionHitRegionEventData region_hit_region_event_data;
        RegionHitBoundaryEventData region_hit_boundary_event_data;
        BlossomShatterEventData blossom_shatter_event_data;
    };
    /// Indicates the type of notification being sent to the mwpm algorithm.
    MwpmEventType event_type;

    MwpmEvent();
    MwpmEvent(RegionHitRegionEventData region_hit_region_event_data);      // NOLINT(google-explicit-constructor)
    MwpmEvent(RegionHitBoundaryEventData region_hit_boundary_event_data);  // NOLINT(google-explicit-constructor)
    MwpmEvent(BlossomShatterEventData blossom_shatter_event_data);         // NOLINT(google-explicit-constructor)
    inline static MwpmEvent no_event() {
        return {};
    }

    bool operator==(const MwpmEvent &rhs) const;
    bool operator!=(const MwpmEvent &rhs) const;
    std::string str() const;
};
std::ostream &operator<<(std::ostream &out, const MwpmEvent &ev);

}  // namespace pm

#endif  // PYMATCHING_FILL_MATCH_MWPM_EVENT_H
