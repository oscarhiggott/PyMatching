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

#ifndef PYMATCHING_FILL_MATCH_QUEUED_EVENT_TRACKER_H
#define PYMATCHING_FILL_MATCH_QUEUED_EVENT_TRACKER_H

#include "pymatching/sparse_blossom/tracker/radix_heap_queue.h"

namespace pm {

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
    inline void set_desired_event(FloodCheckEvent ev, radix_heap_queue<use_validation> &queue) {
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
    /// if doing so was deferred while the earlier event was in the queue.
    template <bool use_validation>
    inline bool dequeue_decision(FloodCheckEvent ev, radix_heap_queue<use_validation> &queue) {
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

#endif  // PYMATCHING_FILL_MATCH_QUEUED_EVENT_TRACKER_H
