#ifndef PYMATCHING_INTS_H
#define PYMATCHING_INTS_H

#include <cstddef>

#include "pymatching/cyclic.h"

namespace pm {

/// This type is used to store observable masks. An observable mask is a bit packed value where the
/// bit 1<<K is set IFF the observable with index K is flipped.
typedef uint64_t obs_int;

/// This type is used to store the weight of an edge.
typedef uint16_t weight_int;

/// This type is used to represent absolute times, accumulated times, and accumulated distances.
/// It is important that it be signed because, for example, it's possible to compute potential
/// collision times that are in the past while considering whether a collision will occur in the
/// future or not.
typedef int32_t cumulative_time_int;

/// This type is used to represent timestamps that are near the current time. To widen a cyclic
/// time to an absolute time, find the absolute time which is closest to the current time and
/// equal to the cyclic time when truncated..
typedef cyclic<uint16_t> cyclic_time_int;

}  // namespace pm

#endif  // PYMATCHING_INTS_H
