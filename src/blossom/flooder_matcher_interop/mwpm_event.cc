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

#include "pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.h"

using namespace pm;

std::ostream &pm::operator<<(std::ostream &out, const RegionHitRegionEventData &dat) {
    out << "{.region1=" << dat.region1;
    out << ", .region2=" << dat.region2;
    out << ", .edge=" << dat.edge;
    out << "}";
    return out;
}

std::ostream &pm::operator<<(std::ostream &out, const RegionHitBoundaryEventData &ev) {
    out << "{.region=" << ev.region;
    out << ", .edge=" << ev.edge;
    out << "}";
    return out;
}

std::ostream &pm::operator<<(std::ostream &out, const BlossomShatterEventData &dat) {
    out << "{.blossom=" << dat.blossom_region;
    out << ", .in_parent_region=" << dat.in_parent_region;
    out << ", .in_child_region=" << dat.in_child_region;
    out << "}";
    return out;
}

std::ostream &pm::operator<<(std::ostream &out, const MwpmEvent &ev) {
    out << "MwpmEvent{.type=";
    switch (ev.event_type) {
        case NO_EVENT:
            out << "NO_EVENT";
            break;
        case REGION_HIT_REGION:
            out << "REGION_HIT_REGION, .dat=" << ev.region_hit_region_event_data;
            break;
        case BLOSSOM_SHATTER:
            out << "BLOSSOM_SHATTER, .dat=" << ev.blossom_shatter_event_data;
            break;
        case REGION_HIT_BOUNDARY:
            out << "REGION_HIT_BOUNDARY, .dat=" << ev.region_hit_boundary_event_data;
            break;
        default:
            throw std::invalid_argument("Unrecognized event type");
    }
    out << "}";
    return out;
}

std::string RegionHitRegionEventData::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

std::string BlossomShatterEventData::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

std::string RegionHitBoundaryEventData::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

std::string MwpmEvent::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

bool RegionHitRegionEventData::operator==(const RegionHitRegionEventData &rhs) const {
    if (region1 == rhs.region1 && region2 == rhs.region2 && edge == rhs.edge) {
        return true;
    }

    // Also equal if reversed.
    return region2 == rhs.region1 && region1 == rhs.region2 && edge == rhs.edge.reversed();
}

bool RegionHitRegionEventData::operator!=(const RegionHitRegionEventData &rhs) const {
    return !(rhs == *this);
}

bool RegionHitBoundaryEventData::operator==(const RegionHitBoundaryEventData &rhs) const {
    return region == rhs.region && edge == rhs.edge;
}

bool RegionHitBoundaryEventData::operator!=(const RegionHitBoundaryEventData &rhs) const {
    return !(rhs == *this);
}

bool BlossomShatterEventData::operator==(const BlossomShatterEventData &rhs) const {
    return blossom_region == rhs.blossom_region && in_parent_region == rhs.in_parent_region &&
           in_child_region == rhs.in_child_region;
}

bool BlossomShatterEventData::operator!=(const BlossomShatterEventData &rhs) const {
    return !(rhs == *this);
}

MwpmEvent::MwpmEvent(RegionHitRegionEventData region_hit_region_event_data)
    : region_hit_region_event_data(region_hit_region_event_data), event_type(REGION_HIT_REGION) {
}
MwpmEvent::MwpmEvent(RegionHitBoundaryEventData region_hit_region_event_data)
    : region_hit_boundary_event_data(region_hit_region_event_data), event_type(REGION_HIT_BOUNDARY) {
}
MwpmEvent::MwpmEvent(BlossomShatterEventData region_hit_region_event_data)
    : blossom_shatter_event_data(region_hit_region_event_data), event_type(BLOSSOM_SHATTER) {
}
MwpmEvent::MwpmEvent() : event_type(NO_EVENT) {
}

bool MwpmEvent::operator==(const MwpmEvent &rhs) const {
    if (event_type != rhs.event_type)
        return false;
    switch (event_type) {
        case NO_EVENT:
            return true;
        case REGION_HIT_REGION:
            return region_hit_region_event_data == rhs.region_hit_region_event_data;
        case REGION_HIT_BOUNDARY:
            return region_hit_boundary_event_data == rhs.region_hit_boundary_event_data;
        case BLOSSOM_SHATTER:
            return blossom_shatter_event_data == rhs.blossom_shatter_event_data;
        default:
            throw std::invalid_argument("Unrecognized event type");
    }
}

bool MwpmEvent::operator!=(const MwpmEvent &rhs) const {
    return !(rhs == *this);
}
