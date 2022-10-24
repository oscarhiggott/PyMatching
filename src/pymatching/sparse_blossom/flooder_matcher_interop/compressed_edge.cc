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

#include "pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.h"

bool pm::CompressedEdge::operator==(const pm::CompressedEdge &rhs) const {
    return loc_from == rhs.loc_from && loc_to == rhs.loc_to && obs_mask == rhs.obs_mask;
}

bool pm::CompressedEdge::operator!=(const pm::CompressedEdge &rhs) const {
    return !(rhs == *this);
}

std::string pm::CompressedEdge::str() const {
    std::stringstream out;
    out << *this;
    return out.str();
}

std::ostream &pm::operator<<(std::ostream &out, const pm::CompressedEdge &edge) {
    out << "CompressedEdge{.obs_mask=" << edge.obs_mask;
    out << ", .loc_from=" << edge.loc_from;
    out << ", .loc_to=" << edge.loc_to;
    out << "}";
    return out;
}
