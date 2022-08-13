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

    bool empty() const {
        return q.empty();
    }

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

    TentativeEvent force_pop() {
        TentativeEvent out{};
        bool b = try_pop(&out);
        if (!b) {
            throw std::invalid_argument("force_pop failed");
        }
        return out;
    }
};

template <>
struct bucket_queue<1> {
    static constexpr size_t BRANCH_FACTOR = 64;
    std::vector<TentativeEvent> heap;
    pm::time_int cur_time;

    explicit bucket_queue(size_t num_buckets) : heap(), cur_time(0) {}

    size_t size() const {
        return heap.size();
    }

    bool satisfies_heap_invariant() const {
        for (size_t k = 0; k < heap.size(); k++) {
            for (size_t c = k*BRANCH_FACTOR + 1; c < std::min(heap.size(), k*BRANCH_FACTOR + BRANCH_FACTOR); c++) {
                if (heap[k].time > heap[c].time) {
                    return false;
                }
            }
        }
        return true;
    }

    void bubble_up(TentativeEvent event, size_t index) {
        auto time = event.time;
        while (index != 0) {
            size_t parent = (index - 1) / BRANCH_FACTOR;
            auto parent_time = heap[parent].time;
            if (parent_time <= time) {
                break;
            }
            heap[index] = heap[parent];
            index = parent;
        }
        heap[index] = event;
    }

    void enqueue(TentativeEvent event) {
        size_t index = heap.size();
        heap.push_back({});
        bubble_up(event, index);
    }

    bool try_pop_once(TentativeEvent *out) {
        if (heap.empty()) {
            return false;
        }
        *out = heap[0];

        size_t index_to_fill = 0;
        while (true) {
            size_t child_start = index_to_fill * BRANCH_FACTOR + 1;
            size_t child_end = std::min(child_start + BRANCH_FACTOR, heap.size());
            if (child_start >= child_end) {
                break;
            }

            size_t best_child = child_start;
            auto best_time = heap[child_start].time;
            for (size_t c = child_start + 1; c < child_end; c++) {
                auto t = heap[c].time;
                if (t < best_time) {
                    best_child = c;
                    best_time = t;
                }
            }

            heap[index_to_fill] = heap[best_child];
            index_to_fill = best_child;
        }

        if (index_to_fill != heap.size() - 1) {
            bubble_up(heap.back(), index_to_fill);
        }
        heap.pop_back();

        return true;
    }

    bool empty() const {
        return heap.empty();
    }

    bool try_pop(TentativeEvent *out) {
        while (try_pop_once(out)) {
            cur_time = out->time;
            if (!out->is_still_valid()) {
                continue;
            }
            return true;
        }
        return false;
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
