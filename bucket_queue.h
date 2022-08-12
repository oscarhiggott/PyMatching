#ifndef PYMATCHING2_BUCKET_QUEUE_H
#define PYMATCHING2_BUCKET_QUEUE_H

#include<vector>
#include "events.h"

namespace pm {

    template <pm::weight_int NUM_BUCKETS>
    class bucket_queue {
    public:
        std::vector<TentativeEvent*> bucket_queue[NUM_BUCKETS];
    };

}

#endif //PYMATCHING2_BUCKET_QUEUE_H
