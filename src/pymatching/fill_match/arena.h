#ifndef PYMATCHING_FILL_MATCH_ARENA_H
#define PYMATCHING_FILL_MATCH_ARENA_H

#include <vector>
#include <iostream>

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
                T* p = (T*)malloc(sizeof(T));
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
            for (T *v : not_in_use) {
                new (v) T();
            }
            for (T *v : to_free) {
                v->~T();
                free(v);
            }
        }
    };
}

#endif