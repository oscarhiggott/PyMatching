#ifndef PYMATCHING2_BUCKET_QUEUE_H
#define PYMATCHING2_BUCKET_QUEUE_H

#include <array>
#include <bit>
#include <iostream>
#include <queue>
#include <vector>

#include "pymatching/cyclic.h"
#include "pymatching/events.h"
#include "pymatching/ints.h"

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
struct bit_bucket_queue {
    std::array<std::vector<TentativeEvent>, sizeof(pm::cyclic_time_int)*8 + 2> bit_buckets;
    pm::cumulative_time_int cur_time;
    size_t _num_enqueued;

    bit_bucket_queue() : cur_time{0}, _num_enqueued(0) {
        // Artificial event just to stop the bucket search.
        bit_buckets.back().push_back(TentativeEvent(cyclic_time_int{0xDEAD}));
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
    void enqueue(TentativeEvent event) {
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
    TentativeEvent dequeue() {
        if (bit_buckets[0].empty()) {
            // Need to refill bucket 0, so we can dequeue from it.

            // Find first non-empty bucket. It has the soonest event.
            size_t b = 1;
            while (bit_buckets[b].empty()) {
                b++;
            }
            if (b == bit_buckets.size() - 1) {
                // We found the fake tail bucket. All real buckets are empty. The queue is empty.
                return TentativeEvent(cyclic_time_int{0});
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
        TentativeEvent result = bit_buckets[0].back();
        bit_buckets[0].pop_back();
        return result;
    }

    /// Lists the sorted events in the queue.
    ///
    /// This method mostly exacts to facilitate testing. It doesn't really make sense to use it
    /// during normal operation.
    std::vector<TentativeEvent> to_vector() const {
        std::vector<TentativeEvent> result;
        for (size_t b = 0; b < bit_buckets.size() - 1; b++) {
            result.insert(result.begin(), bit_buckets[b].begin(), bit_buckets[b].end());
        }
        std::sort(result.begin(),
                  result.end(),
                  [](const TentativeEvent &e1, const TentativeEvent &e2) {
                      return e1.time < e2.time;
                  });
        return result;
    }

    std::string str() const;
};

template <bool use_validation>
std::ostream &operator<<(std::ostream &out, bit_bucket_queue<use_validation> q) {
    out << "bit_bucket_queue {\n";
    out << "    cur_time=" << q.cur_time << "\n";
    for (size_t b = 0; b < q.bit_buckets.size() - 1; b++) {
        auto copy = q.bit_buckets[b];
        if (!copy.empty()) {
            out << "    bucket[" << b << "] {\n";
            std::sort(copy.begin(),
                      copy.end(),
                      [](const TentativeEvent &e1, const TentativeEvent &e2) {
                          return e1.time < e2.time;
                      });
            for (auto &e: copy) {
                out << "        " << e << ",\n";
            }
            out << "    }\n";
        }
    }
    out << "}";
    return out;
}

template <bool use_validation>
std::string bit_bucket_queue<use_validation>::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

/// This class is responsible for ensuring that a "look at me!" event is in the event queue.
///
/// This object also attempts to avoid spamming the event queue with redundant events. If there are
/// two reasons to look at an object at a specific time, it will only put one event into the queue.
/// If an object is already going to be looked at at time T, and it wants to be looked at at time
/// T+2, the look-event for time T+2 will not be enqueued right away. Instead, as part of processing
/// the time T event, the tracker will take care of enqueueing the time T+2 event.
struct QueuedEventTracker {
    cyclic_time_int desired_time{0};
    cyclic_time_int queued_time{0};
    bool has_desired_time{false};
    bool has_queued_time{false};

    /// Resets the tracker to its initial idle state.
    void clear() {
        desired_time = cyclic_time_int{0};
        queued_time = cyclic_time_int{0};
        has_desired_time = {false};
        has_queued_time = {false};
    }

    /// Tells the tracker a desired look-at-me event. The tracker will handle inserting an event
    /// into the event queue, if necessary.
    template <bool use_validation>
    inline void set_desired_event(TentativeEvent ev, bit_bucket_queue<use_validation> &queue) {
        has_desired_time = true;
        desired_time = ev.time;
        if (!has_queued_time || queued_time > ev.time) {
            queued_time = ev.time;
            has_queued_time = true;
            queue.enqueue(ev);
        }
    }

    /// Indicates that it's no longer necessary to look at the object at a later time.
    inline void set_no_desired_event() {
        has_desired_time = false;
    }

    /// Notifies the tracker that a relevant look-at-me event has been dequeued from the event
    /// queue. The result of this method is whether or not to discard the event (due to it being
    /// no longer desired) instead of continuing processing it. Returning false means discard,
    /// returning true means keep. This method also handles requeueing another look-at-me event
    /// if doing os was defered while the earlier event was in the queue.
    template <bool use_validation>
    inline bool dequeue_decision(TentativeEvent ev, bit_bucket_queue<use_validation> &queue) {
        // Only the most recent event this tracker put into the queue is valid. Older events
        // are forgotten and must not be processed, because otherwise an event storm can be
        // created as stale events trigger redundant enqueues.
        if (!has_queued_time || ev.time != queued_time) {
            return false;
        }
        has_queued_time = false;

        // If the event isn't for the CURRENTLY desired look-at-me time, discard it.
        if (!has_desired_time) {
            return false;
        }
        if (ev.time != desired_time) {
            // Requeue the event if the desired time is a little later.
            has_queued_time = true;
            queued_time = desired_time;
            ev.time = desired_time;
            queue.enqueue(ev);
            return false;
        }

        // All systems go! Process away!
        has_desired_time = false;
        return true;
    }
};

}  // namespace pm

#endif  // PYMATCHING2_BUCKET_QUEUE_H
