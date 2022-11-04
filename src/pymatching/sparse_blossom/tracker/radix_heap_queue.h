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

#ifndef PYMATCHING2_BUCKET_QUEUE_H
#define PYMATCHING2_BUCKET_QUEUE_H

#include <algorithm>
#include <array>
#include <bit>
#include <iostream>
#include <queue>
#include <vector>

#include "pymatching/sparse_blossom/ints.h"
#include "pymatching/sparse_blossom/tracker/cyclic.h"
#include "pymatching/sparse_blossom/tracker/flood_check_event.h"

namespace pm {

/// A monotonic priority queue for TentativeEvents.
///
/// The priority queue assumes that times increase monotonically. The caller must not enqueue a
/// time that is cycle-before the time of the last popped event. Time t1 is "cycle-before" time t2
/// iff it takes more increments to get from t1 to t2 than it takes to get from t2 to t1.
///
/// The priority queue works by categorizing upcoming events into groups based on the most
/// significant bit that differs between the event's time and the current time. For example, here
/// are the times that are stored in each bucket assuming the current time is 17:
///
///     (current time = 17)
///     bucket 0: event time in [17]
///     bucket 1: event time in [] (empty because 17 & (1 << 0) is non-zero)
///     bucket 2: event time in [18, 19]
///     bucket 3: event time in [20, 21, 22, 23]
///     bucket 4: event time in [24, 25, 26, 27, 28, 29, 30, 31]
///     bucket 5: event time in [] (empty because 17 & (1 << 4) is non-zero)
///     bucket 6: event time in [32, ..., 63]
///     bucket 7: event time in [64, ..., 127]
///     bucket 8: event time in [128, ..., 255]
///     bucket 9: event time in [256, ..., 511]
///     ...
///     bucket K: event time in [1<<(K-1), (1<<K)-1]
///     ...
///     bucket 31: event time in [1<<30, (1<<31)-1]
///     bucket 32: [invalid for an event to go into this bucket]
///
/// Dequeueing always happens from bucket 0. Whenever bucket 0 is empty, it's refilled by
/// redistributing events from the smallest non-empty bucket (after advancing the time to the
/// minimum time of any event in that bucket).
///
/// The experience of a single event going through the queue is, roughly speaking:
/// - Initially end up in one of the higher buckets.
/// - Sit there for awhile as the lower buckets are drained.
/// - Get redistributed to a slightly lower bucket.
/// - Repeatedly wait and get redistributed to lower and lower buckets until in bucket 0.
/// - Get dequeued out of bucket 0 and yielded as a result.
template <bool use_validation>
struct radix_heap_queue {
    std::array<std::vector<FloodCheckEvent>, sizeof(pm::cyclic_time_int) * 8 + 1> bit_buckets;
    pm::cumulative_time_int cur_time;
    size_t _num_enqueued;

    radix_heap_queue() : cur_time{0}, _num_enqueued(0) {
    }

    size_t size() const {
        return _num_enqueued;
    }

    bool empty() const {
        return _num_enqueued == 0;
    }

    /// Determines which bucket an event with the given time should go into.
    inline size_t cur_bit_bucket_for(cyclic_time_int time) const {
        return std::bit_width((uint64_t)(time.value ^ cyclic_time_int{cur_time}.value));
    }

    /// Adds an event to the priority queue.
    ///
    /// The event MUST NOT be cycle-before the current time.
    void enqueue(FloodCheckEvent event) {
        if (use_validation) {
            if (event.time < cyclic_time_int{cur_time}) {
                std::stringstream ss;
                ss << "Attempted to schedule an event cycle-before the present.\n";
                ss << "    current time: " << cur_time << "\n";
                ss << "    tentative event: " << event << "\n";
                throw std::invalid_argument(ss.str());
            }
        }
        auto target_bucket = cur_bit_bucket_for(event.time);
        bit_buckets[target_bucket].push_back(event);
        _num_enqueued++;
    }

    /// Checks if all events are in the correct bucket.
    bool satisfies_invariants() const {
        for (size_t b = 0; b < bit_buckets.size() - 1; b++) {
            for (size_t k = 0; k < bit_buckets[b].size(); k++) {
                if (cur_bit_bucket_for(bit_buckets[b][k].time) != b) {
                    return false;
                }
            }
        }
        return true;
    }

    /// Dequeues the next event.
    ///
    /// If the queue is empty, a tentative event with type NO_TENTATIVE_EVENT is returned.
    FloodCheckEvent dequeue() {
        if (_num_enqueued == 0)
            return FloodCheckEvent(cyclic_time_int{0});
        if (bit_buckets[0].empty()) {
            // Need to refill bucket 0, so we can dequeue from it.

            // Find first non-empty bucket. It has the soonest event.
            size_t b = 1;
            while (bit_buckets[b].empty()) {
                b++;
            }

            if (b == 1) {
                // Special case:
                std::swap(bit_buckets[0], bit_buckets[1]);
                cur_time++;
            } else {
                auto &source_bucket = bit_buckets[b];

                // Advance time to the minimum time in the bucket.
                decltype(cyclic_time_int::value) min_time = source_bucket[0].time.value;
                for (size_t k = 1; k < source_bucket.size(); k++) {
                    min_time = std::min(min_time, source_bucket[k].time.value);
                }
                cur_time = cyclic_time_int{min_time}.widen_from_nearby_reference(cur_time);

                // Redistribute the contents of the bucket.
                for (auto &e : source_bucket) {
                    auto target_bucket = cur_bit_bucket_for(e.time);
                    bit_buckets[target_bucket].push_back(e);
                }
                source_bucket.clear();
            }
        }

        _num_enqueued--;
        FloodCheckEvent result = bit_buckets[0].back();
        bit_buckets[0].pop_back();
        return result;
    }

    /// Lists the sorted events in the queue.
    ///
    /// This method mostly exacts to facilitate testing. It doesn't really make sense to use it
    /// during normal operation.
    std::vector<FloodCheckEvent> to_vector() const {
        std::vector<FloodCheckEvent> result;
        for (size_t b = 0; b < bit_buckets.size() - 1; b++) {
            result.insert(result.begin(), bit_buckets[b].begin(), bit_buckets[b].end());
        }
        std::sort(result.begin(), result.end(), [](const FloodCheckEvent &e1, const FloodCheckEvent &e2) {
            return e1.time < e2.time;
        });
        return result;
    }

    std::string str() const;

    /// Clear all remaining events from the queue
    void clear();

    /// Clear the queue and set cur_time = 0
    void reset();
};

template <bool use_validation>
std::ostream &operator<<(std::ostream &out, radix_heap_queue<use_validation> q) {
    out << "bit_bucket_queue {\n";
    out << "    cur_time=" << q.cur_time << "\n";
    for (size_t b = 0; b < q.bit_buckets.size() - 1; b++) {
        auto copy = q.bit_buckets[b];
        if (!copy.empty()) {
            out << "    bucket[" << b << "] {\n";
            std::sort(copy.begin(), copy.end(), [](const FloodCheckEvent &e1, const FloodCheckEvent &e2) {
                return e1.time < e2.time;
            });
            for (auto &e : copy) {
                out << "        " << e << ",\n";
            }
            out << "    }\n";
        }
    }
    out << "}";
    return out;
}

template <bool use_validation>
std::string radix_heap_queue<use_validation>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

template <bool use_validation>
void radix_heap_queue<use_validation>::clear() {
    size_t b = 0;
    while (!empty()) {
        _num_enqueued -= bit_buckets[b].size();
        bit_buckets[b].clear();
        b++;
    }
}

template <bool use_validation>
void radix_heap_queue<use_validation>::reset() {
    clear();
    cur_time = 0;
}

}  // namespace pm

#endif  // PYMATCHING2_BUCKET_QUEUE_H
