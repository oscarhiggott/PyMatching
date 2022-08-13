#ifndef PYMATCHING2_BUCKET_QUEUE_H
#define PYMATCHING2_BUCKET_QUEUE_H

#include<vector>
#include "events.h"

namespace pm {

    template <pm::weight_int NUM_QUEUE_BUCKETS>
    class bucket_queue {
    public:
        std::vector<TentativeEvent> buckets_array[NUM_QUEUE_BUCKETS];
        pm::weight_int current_bucket;
        bucket_queue();
        template< class... Args >
        void emplace( Args&&... args );
        void push(const TentativeEvent& event);
        const TentativeEvent& top() const;
        void pop();
        bool empty() const;
    };

    template<pm::weight_int NUM_QUEUE_BUCKETS>
    bucket_queue<NUM_QUEUE_BUCKETS>::bucket_queue() : current_bucket(0) {}

    template<pm::weight_int NUM_QUEUE_BUCKETS>
    void bucket_queue<NUM_QUEUE_BUCKETS>::push(const TentativeEvent &event) {
        pm::weight_int bucket = event.time % NUM_QUEUE_BUCKETS;
        buckets_array[bucket].push_back(event);
    }

}

#endif //PYMATCHING2_BUCKET_QUEUE_H
