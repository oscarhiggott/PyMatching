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

#include "pymatching/sparse_blossom/tracker/cyclic.h"

#include <gtest/gtest.h>

using namespace pm;

TEST(cyclic, comparison) {
    for (size_t offset = 0; offset < 256; offset++) {
        cyclic<uint8_t> c0{(uint8_t)(0 + offset)};
        cyclic<uint8_t> c1{(uint8_t)(1 + offset)};
        cyclic<uint8_t> c127{(uint8_t)(127 + offset)};
        cyclic<uint8_t> c128{(uint8_t)(128 + offset)};
        cyclic<uint8_t> c129{(uint8_t)(129 + offset)};
        cyclic<uint8_t> c255{(uint8_t)(255 + offset)};

        ASSERT_TRUE(c0 >= c0);
        ASSERT_FALSE(c0 > c0);
        ASSERT_FALSE(c0 < c0);
        ASSERT_TRUE(c0 <= c0);
        ASSERT_TRUE(c0 == c0);
        ASSERT_FALSE(c0 != c0);

        ASSERT_FALSE(c0 >= c1);
        ASSERT_FALSE(c0 > c1);
        ASSERT_TRUE(c0 < c1);
        ASSERT_TRUE(c0 <= c1);
        ASSERT_FALSE(c0 == c1);
        ASSERT_TRUE(c0 != c1);

        ASSERT_FALSE(c0 >= c127);
        ASSERT_FALSE(c0 > c127);
        ASSERT_TRUE(c0 < c127);
        ASSERT_TRUE(c0 <= c127);
        ASSERT_FALSE(c0 == c127);
        ASSERT_TRUE(c0 != c127);

        ASSERT_FALSE(c0 >= c128);
        ASSERT_FALSE(c0 > c128);
        ASSERT_FALSE(c0 < c128);
        ASSERT_FALSE(c0 <= c128);
        ASSERT_FALSE(c0 == c128);
        ASSERT_TRUE(c0 != c128);

        ASSERT_TRUE(c0 >= c129);
        ASSERT_TRUE(c0 > c129);
        ASSERT_FALSE(c0 < c129);
        ASSERT_FALSE(c0 <= c129);
        ASSERT_FALSE(c0 == c129);
        ASSERT_TRUE(c0 != c129);

        ASSERT_TRUE(c0 >= c255);
        ASSERT_TRUE(c0 > c255);
        ASSERT_FALSE(c0 < c255);
        ASSERT_FALSE(c0 <= c255);
        ASSERT_FALSE(c0 == c255);
        ASSERT_TRUE(c0 != c255);

        ASSERT_TRUE(c1 >= c0);
        ASSERT_TRUE(c1 > c0);
        ASSERT_FALSE(c1 < c0);
        ASSERT_FALSE(c1 <= c0);
        ASSERT_FALSE(c1 == c0);
        ASSERT_TRUE(c1 != c0);

        ASSERT_TRUE(c127 >= c0);
        ASSERT_TRUE(c127 > c0);
        ASSERT_FALSE(c127 < c0);
        ASSERT_FALSE(c127 <= c0);
        ASSERT_FALSE(c127 == c0);
        ASSERT_TRUE(c127 != c0);

        ASSERT_FALSE(c128 >= c0);
        ASSERT_FALSE(c128 > c0);
        ASSERT_FALSE(c128 < c0);
        ASSERT_FALSE(c128 <= c0);
        ASSERT_FALSE(c128 == c0);
        ASSERT_TRUE(c128 != c0);

        ASSERT_FALSE(c129 >= c0);
        ASSERT_FALSE(c129 > c0);
        ASSERT_TRUE(c129 < c0);
        ASSERT_TRUE(c129 <= c0);
        ASSERT_FALSE(c129 == c0);
        ASSERT_TRUE(c129 != c0);

        ASSERT_FALSE(c255 >= c0);
        ASSERT_FALSE(c255 > c0);
        ASSERT_TRUE(c255 < c0);
        ASSERT_TRUE(c255 <= c0);
        ASSERT_FALSE(c255 == c0);
        ASSERT_TRUE(c255 != c0);
    }
}

TEST(cyclic, arithmetic) {
    cyclic<uint16_t> c{0};
    c += 1;
    ASSERT_EQ(c, 1);
    ASSERT_EQ(c++, 1);
    ASSERT_EQ(c, 2);
    ASSERT_EQ(++c, 3);
    ASSERT_EQ(c, 3);
    ASSERT_EQ(c--, 3);
    ASSERT_EQ(c, 2);
    ASSERT_EQ(--c, 1);
    ASSERT_EQ(c, 1);
    c -= 5;
    ASSERT_EQ(c, (uint16_t)-4);

    cyclic<uint16_t> d{8};
    ASSERT_EQ(c + d, 4);
    ASSERT_EQ(d + c, 4);
    ASSERT_EQ(c - d, (uint16_t)-12);
    ASSERT_EQ(d - c, 12);
}

TEST(cyclic, widen_from_nearby_reference) {
    cyclic<uint8_t> v{130};
    ASSERT_EQ(v.widen_from_nearby_reference(2048u + 130u), 130u + 2048u);
    ASSERT_EQ(v.widen_from_nearby_reference(2048u + 100u), 130u + 2048u);
    ASSERT_EQ(v.widen_from_nearby_reference(2048u + 130u - 127u), 130u + 2048u);
    ASSERT_EQ(v.widen_from_nearby_reference(2048u + 130u + 127u), 130u + 2048u);
    ASSERT_EQ(v.widen_from_nearby_reference(2048u + 130u - 129u), 130u + 2048u - 256u);
    ASSERT_EQ(v.widen_from_nearby_reference(2048u + 130u + 129u), 130u + 2048u + 256u);

    ASSERT_EQ(v.widen_from_nearby_reference(-2048 + 130), 130 - 2048);
    ASSERT_EQ(v.widen_from_nearby_reference(-2048 + 100), 130 - 2048);
    ASSERT_EQ(v.widen_from_nearby_reference(-2048 + 130 - 127), 130 - 2048);
    ASSERT_EQ(v.widen_from_nearby_reference(-2048 + 130 + 127), 130 - 2048);
    ASSERT_EQ(v.widen_from_nearby_reference(-2048 + 130 - 129), 130 - 2048 - 256u);
    ASSERT_EQ(v.widen_from_nearby_reference(-2048 + 130 + 129), 130 - 2048 + 256u);

    ASSERT_EQ(cyclic<uint16_t>{16}.widen_from_nearby_reference(int32_t{65535}), 65536 + 16);
    ASSERT_EQ(cyclic<uint16_t>{16}.widen_from_nearby_reference(int32_t{-65535}), -65536 + 16);
}
