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

#ifndef PYMATCHING_FILL_MATCH_ARENA_H
#define PYMATCHING_FILL_MATCH_ARENA_H

#include <algorithm>
#include <iostream>
#include <vector>

namespace pm {

/// World's simplest bulk memory owner.
///
/// Memory allocated by the arena is free'd when the arena is destructed.
template <typename T>
struct Arena {
    std::vector<T *> allocated;
    std::vector<T *> available;

    Arena() : allocated(), available() {
    }
    Arena(const Arena &) = delete;
    Arena(Arena &&other) : allocated(std::move(other.allocated)), available(std::move(other.available)) {
    }

    T *alloc_unconstructed() {
        if (available.empty()) {
            T *p = (T *)malloc(sizeof(T));
            allocated.push_back(p);
            available.push_back(p);
        }
        T *result = available.back();
        available.pop_back();
        return result;
    }

    T *alloc_default_constructed() {
        T *result = alloc_unconstructed();
        new (result) T();
        return result;
    }

    void del(T *p) {
        available.push_back(p);
        p->~T();
    }

    ~Arena() {
        std::vector<T *> to_free = std::move(allocated);
        std::vector<T *> not_in_use = std::move(available);

        // Destruct the objects that were still in use.
        std::sort(to_free.begin(), to_free.end());
        std::sort(not_in_use.begin(), not_in_use.end());
        size_t kf = 0;
        size_t kn = 0;
        while (kf < to_free.size()) {
            if (kn < not_in_use.size() && not_in_use[kn] == to_free[kf]) {
                kf++;
                kn++;
            } else {
                to_free[kf]->~T();
                kf++;
            }
        }

        // Free all allocated memory.
        for (T *v : to_free) {
            free(v);
        }
    }
};
}  // namespace pm

#endif