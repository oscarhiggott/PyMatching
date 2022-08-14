#ifndef PYMATCHING2_BUCKET_QUEUE_H
#define PYMATCHING2_BUCKET_QUEUE_H

#include <vector>
#include <queue>
#include <bit>
#include <iostream>

#include "pymatching/events.h"

namespace pm {

/// A priority queue for TentativeEvents.
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
struct bit_bucket_queue {
    std::array<std::vector<TentativeEvent>, sizeof(pm::time_int)*8 + 2> bit_buckets;
    pm::time_int cur_time;
    size_t _num_enqueued;

    bit_bucket_queue() : cur_time(0), _num_enqueued(0) {
        // Artificial event just to stop the bucket search.
        bit_buckets.back().push_back(TentativeEvent(-1, 0xDEAD));
    }

    size_t size() const {
        return _num_enqueued;
    }

    bool empty() const {
        return _num_enqueued == 0;
    }

    /// Determines which bucket an event with the given time should go into.
    inline size_t cur_bit_bucket_for(pm::time_int time) const {
        return std::bit_width((uint32_t)(time ^ cur_time));
    }

    /// Adds an event to the priority queue.
    ///
    /// The event MUST NOT be cycle-before the current time.
    void enqueue(TentativeEvent event) {
        if (use_validation) {
            uint32_t d = (uint32_t) event.time - (uint32_t) cur_time;
            if (d >= 1 << 31) {
                std::stringstream ss;
                if (event.time < cur_time) {
                    ss << "Attempted to schedule an event cycle-before the present.\n";
                    ss << "    current time: " << cur_time << "\n";
                    ss << "    tentative event: " << event << "\n";
                }
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

    /// Dequeues the next event, regardless of whether it is invalidated.
    bool try_pop_any(TentativeEvent *out) {
        if (bit_buckets[0].empty()) {
            // Need to refill bucket 0, so we can dequeue from it.

            // Find first non-empty bucket. It has the soonest event.
            size_t b = 1;
            while (bit_buckets[b].empty()) {
                b++;
            }
            if (b == bit_buckets.size() - 1) {
                // We found the fake tail bucket. All real buckets are empty. The queue is empty.
                return false;
            }

            if (b == 1) {
                // Special case:
                std::swap(bit_buckets[0], bit_buckets[1]);
                cur_time++;
            } else {
                auto &source_bucket = bit_buckets[b];

                cur_time = source_bucket[0].time;
                for (size_t k = 1; k < source_bucket.size(); k++) {
                    cur_time = std::min(cur_time, source_bucket[k].time);
                }

                for (size_t k = 0; k < source_bucket.size(); k++) {
                    auto target_bucket = cur_bit_bucket_for(source_bucket[k].time);
                    if (target_bucket != b) {
                        bit_buckets[target_bucket].push_back(source_bucket[k]);
                        source_bucket[k] = source_bucket.back();
                        source_bucket.pop_back();
                        k--;
                    }
                }
            }
        }

        _num_enqueued--;
        *out = bit_buckets[0].back();
        bit_buckets[0].pop_back();
        return true;
    }

    /// Dequeues the next event, skipping over invalidated events.
    ///
    /// Args:
    ///     out: Where the dequeued event is written, unless the queue is empty.
    ///
    /// Returns:
    ///     true: Event successfully dequeued.
    ///     false: Queue was empty.
    bool try_pop_valid(TentativeEvent *out) {
        while (try_pop_any(out)) {
            cur_time = out->time;
            if (!out->is_still_valid()) {
                continue;
            }
            return true;
        }
        return false;
    }

    /// Dequeues the next event, skipping over invalidated events.
    ///
    /// Returns:
    ///     The dequeued event.
    ///
    /// Raises:
    ///     std::invalid_argument: The queue was empty.
    TentativeEvent force_pop_valid() {
        TentativeEvent out{};
        bool b = try_pop_valid(&out);
        if (!b) {
            throw std::invalid_argument("force_pop_valid failed");
        }
        return out;
    }
};

}  // namespace pm

#endif  // PYMATCHING2_BUCKET_QUEUE_H
