#include "chunk_list.h"

#include <gtest/gtest.h>

TEST(ChunkList, push_anywhere) {
    pm::ChunkList<uint16_t, 4> list;
    ASSERT_EQ(list.storage.num_items, 0);
    ASSERT_EQ(list.storage.items, (std::array<uint16_t, 4>{0, 0, 0, 0}));
    ASSERT_FALSE((bool)list.storage.next);
    ASSERT_EQ(list.to_vector(), (std::vector<uint16_t>{}));

    list.push_anywhere(5);
    ASSERT_EQ(list.storage.num_items, 1);
    ASSERT_EQ(list.storage.items, (std::array<uint16_t, 4>{5, 0, 0, 0}));
    ASSERT_FALSE((bool)list.storage.next);
    ASSERT_EQ(list.to_vector(), (std::vector<uint16_t>{5}));

    list.push_anywhere(6);
    list.push_anywhere(7);
    list.push_anywhere(8);
    ASSERT_EQ(list.storage.num_items, 4);
    ASSERT_EQ(list.storage.items, (std::array<uint16_t, 4>{5, 6, 7, 8}));
    ASSERT_FALSE((bool)list.storage.next);
    ASSERT_EQ(list.to_vector(), (std::vector<uint16_t>{5, 6, 7, 8}));

    list.push_anywhere(9);
    ASSERT_EQ(list.storage.num_items, 1);
    ASSERT_EQ(list.storage.items, (std::array<uint16_t, 4>{9, 6, 7, 8}));
    ASSERT_TRUE((bool)list.storage.next);
    ASSERT_EQ(list.storage.next->num_items, 4);
    ASSERT_EQ(list.storage.next->items, (std::array<uint16_t, 4>{5, 6, 7, 8}));
    ASSERT_EQ(list.to_vector(), (std::vector<uint16_t>{9, 5, 6, 7, 8}));
    ASSERT_FALSE((bool)list.storage.next->next);
}

TEST(ChunkList, iter_drain_happens_to_reverse_order) {
    pm::ChunkList<uint16_t, 4> list;
    std::set<uint16_t> expected;
    for (uint16_t k = 0; k < 20; k++) {
        list.push_anywhere(k);
    }
    std::vector<uint16_t> drained;
    list.iter_drain_reentrant_safe([&](uint16_t value) {
        drained.push_back(value);
    });
    ASSERT_EQ(
        drained,
        (std::vector<uint16_t>{
            19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
        }));
}

TEST(ChunkList, iter_drain_while_reentering) {
    pm::ChunkList<uint16_t, 4> list;
    list.push_anywhere(0);
    std::vector<uint16_t> drained;
    list.iter_drain_reentrant_safe([&](uint16_t value) {
        if (value < 10) {
            list.push_anywhere(value * 2 + 1);
            list.push_anywhere(value * 2 + 2);
        }
        drained.push_back(value);
    });
    ASSERT_EQ(drained.size(), 21);
    ASSERT_EQ(
        drained, (std::vector<uint16_t>{0, 2, 6, 14, 13, 5, 12, 11, 1, 4, 10, 9, 20, 19, 3, 8, 18, 17, 7, 16, 15}));
}
