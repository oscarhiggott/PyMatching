#ifndef PYMATCHING2_BUCKET_QUEUE_H
#define PYMATCHING2_BUCKET_QUEUE_H

#include <vector>
#include <queue>

#include "pymatching/events.h"
#include "pymatching/chunk_list.h"

namespace pm {

template <size_t CHUNK_SIZE>
struct bucket_queue {
    std::vector<pm::ChunkList<TentativeEvent, CHUNK_SIZE>> buckets;
    size_t _num_enqueued;
    pm::time_int cur_time;

    explicit bucket_queue(size_t num_buckets) : buckets(num_buckets), _num_enqueued(0), cur_time(0) {}

    size_t size() const {
        return _num_enqueued;
    }

    void enqueue(TentativeEvent event) {
//        if (event->time < cur_time || event->time >= cur_time + buckets.size()) {
//            throw std::invalid_argument("event->time is out of bucket range");
//        }
        buckets[event.time % buckets.size()].push_anywhere(event);
        _num_enqueued++;
    }

    bool try_pop(TentativeEvent *out) {
        while (true) {
            auto &bucket = buckets[cur_time % buckets.size()];
            if (bucket.try_pop(out)) {
                _num_enqueued--;
                if (!out->is_still_valid()) {
                    continue;
                }
                return true;
            }
            if (_num_enqueued == 0) {
                return false;
            }
            cur_time++;
        }
    }

    bool empty() const {
        return _num_enqueued == 0;
    }
    TentativeEvent force_pop() {
        TentativeEvent out{};
        bool b = try_pop(&out);
        if (!b) {
            throw std::invalid_argument("force_pop failed");
        }
        return out;
    }
};

class bucket_queue_compare_helper {
   public:
    bool operator()(TentativeEvent a, TentativeEvent b) {
        return a > b;
    }
};

template <>
struct bucket_queue<0> {
    std::priority_queue<TentativeEvent, std::vector<TentativeEvent>, bucket_queue_compare_helper> q;
    pm::time_int cur_time;

    explicit bucket_queue(size_t num_buckets) : q(), cur_time(0) {}

    size_t size() const {
        return q.size();
    }

    void enqueue(TentativeEvent event) {
        q.push(event);
    }

    bool try_pop(TentativeEvent *out) {
        while (!q.empty()) {
            *out = q.top();
            q.pop();
            cur_time = out->time;
            if (!out->is_still_valid()) {
                continue;
            }
            return true;
        }
        return false;
    }

    bool empty() const {
        return q.empty();
    }

    TentativeEvent force_pop() {
        TentativeEvent out{};
        bool b = try_pop(&out);
        if (!b) {
            throw std::invalid_argument("force_pop failed");
        }
        return out;
    }
};

}  // namespace pm

#endif  // PYMATCHING2_BUCKET_QUEUE_H
