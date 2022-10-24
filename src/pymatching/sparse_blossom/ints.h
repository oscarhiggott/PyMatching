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

#ifndef PYMATCHING_INTS_H
#define PYMATCHING_INTS_H

#include <cstddef>

#include "pymatching/sparse_blossom/tracker/cyclic.h"

namespace pm {

/// This type is used to store observable masks. An observable mask is a bit packed value where the
/// bit 1<<K is set IFF the observable with index K is flipped.
typedef uint64_t obs_int;

/// This type is used to store the weight of an edge.
typedef uint32_t weight_int;

/// This type is used to store the potentially-negative weight of an edge.
/// It is used when loading graphs in order to support negative edge weights.
/// However, negative edge weights are handled in pre- and post-processing, rather
/// than in the blossom algorithm itself, so the `weight_int' type is used instead for data
/// structures used within the blossom algorithm.
typedef int32_t signed_weight_int;

/// This type is used to represent absolute times, accumulated times, and accumulated distances.
/// It is important that it be signed because, for example, it's possible to compute potential
/// collision times that are in the past while considering whether a collision will occur in the
/// future or not.
typedef int64_t cumulative_time_int;

/// This type is used to represent the total weight of the MWPM solution. It is important that it
/// is 64-bit since in general we expect the total solution weight to grow linearly with the number of
/// edges in the graph. It is also important that it is signed, since negative weight edges are permitted
/// as input to the decoder (albeit not within the blossom algorithm itself).
typedef int64_t total_weight_int;

/// This type is used to represent timestamps that are near the current time. To widen a cyclic
/// time to an absolute time, find the absolute time which is closest to the current time and
/// equal to the cyclic time when truncated..
typedef cyclic<uint32_t> cyclic_time_int;

}  // namespace pm

#endif  // PYMATCHING_INTS_H
