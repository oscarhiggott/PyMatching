#ifndef PYMATCHING2_FIXED_LENGTH_VECTOR_H
#define PYMATCHING2_FIXED_LENGTH_VECTOR_H

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace pm {

template <typename T, size_t N>
struct ChunkListStorage {
    size_t num_items;
    std::unique_ptr<ChunkListStorage<T, N>> next;
    std::array<T, N> items;
};

template <typename T, size_t N>
struct ChunkList {
    ChunkListStorage<T, N> storage;

    ChunkList() : storage{0} {
    }

    inline void push_anywhere(T item) {
        if (storage.num_items >= N) {
            auto moved = new ChunkListStorage<T, N>(std::move(storage));
            storage.next.reset(moved);
            storage.num_items = 0;
        }
        storage.items[storage.num_items] = item;
        storage.num_items++;
    }

    std::vector<T> to_vector() const {
        std::vector<T> result;
        const ChunkListStorage<T, N> *s = &storage;
        while (s != nullptr) {
            for (size_t k = 0; k < s->num_items; k++) {
                result.push_back(s->items[k]);
            }
            s = s->next.get();
        }
        return result;
    }

    /// Empties the list, passing each item to a callback.
    /// It is safe for the callback to add items into the list, even as it is being drained.
    template <typename CALLBACK>
    inline void iter_drain_reentrant_safe(CALLBACK callback) {
        while (true) {
            if (storage.num_items == 0) {
                if (!storage.next) {
                    break;
                }
                storage.items = storage.next->items;
                storage.num_items = storage.next->num_items;
                std::unique_ptr<ChunkListStorage<T, N>> hop_ptr = std::move(storage.next->next);
                storage.next = std::move(hop_ptr);
            }
            storage.num_items--;
            T item = std::move(storage.items[storage.num_items]);
            // Don't combine reading the item with invoking the callback.
            // Otherwise the callback may take a reference into the array, and the
            // array slot may be overwritten invalidating that reference while the
            // callback continues to use it.
            callback(item);
        }
    }
};
}  // namespace pm

#endif  // PYMATCHING2_FIXED_LENGTH_VECTOR_H
