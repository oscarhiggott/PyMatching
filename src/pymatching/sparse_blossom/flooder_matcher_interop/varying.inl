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

#include "varying.h"

template <typename T>
inline Varying<T>::Varying() = default;

template <typename T>
inline Varying<T>::Varying(T data) : data(data) {
}

template <typename T>
inline T Varying<T>::get_distance_at_time(T time) const {
    if (data & 1) {
        return (data >> 2) + time;
    } else if (data & 2) {
        return (data >> 2) - time;
    } else {
        return data >> 2;
    }
}

template <typename T>
inline T Varying<T>::get_growing_distance_at_time(T time) const {
    return (data >> 2) + time;
}

template <typename T>
inline T Varying<T>::get_shrinking_distance_at_time(T time) const {
    return (data >> 2) - time;
}

template <typename T>
inline T Varying<T>::time_of_x_intercept() const {
    if (data & 1) {
        return (data >> 2) * -1;
    } else if (data & 2) {
        return data >> 2;
    } else {
        throw std::invalid_argument("Varying must be growing or shrinking to have an x intercept.");
    }
}

template <typename T>
T Varying<T>::time_of_x_intercept_for_growing() const {
    // Assumes Varying is growing
    return -(data >> 2);
}

template <typename T>
T Varying<T>::time_of_x_intercept_for_shrinking() const {
    return data >> 2;
}

template <typename T>
inline T Varying<T>::time_of_x_intercept_when_added_to(Varying<T> other) const {
    // Assumes both this and other are growing or frozen. Cannot have both frozen, or any shrinking.
    T time_with_unit_slope = -(data >> 2) - (other.data >> 2);
    if ((this->data & 1) && (other.data & 1)) {
        return time_with_unit_slope >> 1;
    } else {
        // By assumption on the input, one must be growing, and the other frozen.
        return time_with_unit_slope;
    }
}

template <typename T>
inline T Varying<T>::time_of_x_intercept_when_added_to_giving_unit_slope(Varying<T> other) const {
    // Assumes both this and other are growing or frozen. Cannot have both frozen, or any shrinking.
    return -(data >> 2) - (other.data >> 2);
}

template <typename T>
bool Varying<T>::colliding_with(Varying<T> other) const {
    return ((data | other.data) & 3) == 1;
}

template <typename T>
inline Varying<T> Varying<T>::then_growing_at_time(T time_of_change) const {
    T d = get_distance_at_time(time_of_change);
    return pm::Varying<T>(((d - time_of_change) << 2) + 1);
}

template <typename T>
inline Varying<T> Varying<T>::then_shrinking_at_time(T time_of_change) const {
    T d = get_distance_at_time(time_of_change);
    return pm::Varying<T>(((d + time_of_change) << 2) + 2);
}

template <typename T>
inline Varying<T> Varying<T>::then_frozen_at_time(T time_of_change) const {
    T d = get_distance_at_time(time_of_change);
    return pm::Varying<T>(d << 2);
}

template <typename T>
inline std::string Varying<T>::str() const {
    std::stringstream s;
    s << *this;
    return s.str();
}

template <typename T>
inline Varying<T> Varying<T>::growing_varying_with_zero_distance_at_time(T time) {
    return Varying<T>(((-time) << 2) + 1);
}

template <typename T>
Varying<T> Varying<T>::from_base_and_growth(T base, int8_t growth) {
    Varying<T> result{base << 2};
    if (growth == +1) {
        result.data |= 1;
    } else if (growth == -1) {
        result.data |= 2;
    } else if (growth != 0) {
        throw std::invalid_argument("growth must be 0 or +1 or -1");
    }
    return result;
}

template <typename T>
Varying<T> Varying<T>::growing_value_at_time(T y, T t) {
    return Varying<T>::from_base_and_growth(y - t, +1);
}

template <typename T>
Varying<T> Varying<T>::shrinking_value_at_time(T y, T t) {
    return Varying<T>::from_base_and_growth(y + t, -1);
}

template <typename T>
Varying<T> Varying<T>::frozen(T base) {
    return Varying<T>::from_base_and_growth(base, 0);
}

template <typename T>
inline std::ostream &operator<<(std::ostream &os, Varying<T> varying) {
    os << (varying.data >> 2);
    if (varying.data & 1) {
        os << " + t";
    } else if (varying.data & 2) {
        os << " - t";
    }
    return os;
}

template <typename T>
inline bool Varying<T>::operator==(Varying<T> rhs) const {
    return data == rhs.data;
}

template <typename T>
inline bool Varying<T>::operator!=(Varying<T> rhs) const {
    return data != rhs.data;
}

template <typename T>
inline Varying<T> Varying<T>::operator+(T offset) const {
    return Varying<T>(data + (offset << 2));
}

template <typename T>
inline Varying<T> Varying<T>::operator-(T offset) const {
    return pm::Varying<T>(data - (offset << 2));
}

template <typename T>
inline Varying<T> &Varying<T>::operator+=(T offset) {
    data += offset << 2;
    return *this;
}

template <typename T>
inline Varying<T> &Varying<T>::operator-=(T offset) {
    data -= offset << 2;
    return *this;
}

template <typename T>
T Varying<T>::y_intercept() const {
    return data >> 2;
}

template <typename T>
bool Varying<T>::is_shrinking() const {
    return data & 2;
}

template <typename T>
bool Varying<T>::is_frozen() const {
    return !(data & 3);
}

template <typename T>
bool Varying<T>::is_growing() const {
    return data & 1;
}
