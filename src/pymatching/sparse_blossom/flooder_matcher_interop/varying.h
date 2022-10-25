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

#ifndef PYMATCHING2_VARYING_H
#define PYMATCHING2_VARYING_H

#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include "pymatching/sparse_blossom/ints.h"

namespace pm {

/// A Varying is a value growing linearly with time.
/// The only allowed growth rates are +1, 0, or -1.
template <typename T>
struct Varying {
    /// Bit packed value with top 2 bits storing the slope of the line and
    /// the remaining bits storing the Y intercept of the line.
    T data{};

    Varying();
    explicit Varying(T data);

    T get_distance_at_time(T time) const;
    T get_shrinking_distance_at_time(T time) const;
    T get_growing_distance_at_time(T time) const;
    T time_of_x_intercept() const;
    T time_of_x_intercept_when_added_to(Varying<T> other) const;
    T time_of_x_intercept_when_added_to_giving_unit_slope(Varying<T> other) const;
    T time_of_x_intercept_for_growing() const;
    T time_of_x_intercept_for_shrinking() const;
    bool is_growing() const;
    bool is_shrinking() const;
    bool is_frozen() const;
    bool colliding_with(Varying<T> other) const;
    Varying<T> then_growing_at_time(T time_of_change) const;
    Varying<T> then_shrinking_at_time(T time_of_change) const;
    Varying<T> then_frozen_at_time(T time_of_change) const;
    T y_intercept() const;
    std::string str() const;
    static Varying<T> growing_varying_with_zero_distance_at_time(T time);
    static Varying<T> from_base_and_growth(T base, int8_t growth);
    static Varying<T> frozen(T base);
    static Varying<T> growing_value_at_time(T y, T t);
    static Varying<T> shrinking_value_at_time(T y, T t);
    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, Varying<U> varying);
    bool operator==(Varying<T> rhs) const;
    bool operator!=(Varying<T> rhs) const;
    Varying<T> operator+(T offset) const;
    Varying<T> operator-(T offset) const;
    Varying<T> &operator+=(T offset);
    Varying<T> &operator-=(T offset);
};

using Varying32 = Varying<int32_t>;
using Varying64 = Varying<int64_t>;
using VaryingCT = Varying<cumulative_time_int>;

#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.inl"

}  // namespace pm

#endif  // PYMATCHING2_VARYING_H
