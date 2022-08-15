#ifndef PYMATCHING2_VARYING_H
#define PYMATCHING2_VARYING_H

#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

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
    T time_of_x_intercept();
    T time_of_x_intercept_when_added_to(Varying<T> other);
    T time_of_x_intercept_for_growing();
    T time_of_x_intercept_for_shrinking();
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
    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, Varying<U> varying);
    bool operator==(Varying<T> rhs) const;
    bool operator!=(Varying<T> rhs) const;
    Varying<T> operator+(T offset) const;
    Varying<T> operator-(T offset) const;
    Varying<T> &operator+=(T offset);
    Varying<T> &operator-=(T offset);

    void inplace_freeze_then_add(Varying<T> other);
};

using Varying32 = Varying<int32_t>;
using Varying64 = Varying<int64_t>;

#include "pymatching/fill_match/flooder/varying.inl"

}  // namespace pm

#endif  // PYMATCHING2_VARYING_H
