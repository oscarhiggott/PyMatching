#ifndef PYMATCHING_FILL_MATCH_ARENA_H
#define PYMATCHING_FILL_MATCH_ARENA_H

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