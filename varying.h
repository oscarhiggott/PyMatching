#ifndef PYMATCHING2_VARYING_H
#define PYMATCHING2_VARYING_H

#include <cstdint>
#include <string>
#include <sstream>
#include <stdexcept>
#include <limits>


namespace pm {
    typedef int32_t time_int;

    /// A Varying is a value growing linearly with time.
    /// The only allowed growth rates are +1, 0, or -1.
    template<typename T>
    class Varying {
    public:
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
        bool colliding_with(Varying <T> other) const;
        Varying<T> then_growing_at_time(T time_of_change) const;
        Varying<T> then_shrinking_at_time(T time_of_change) const;
        Varying<T> then_frozen_at_time(T time_of_change) const;
        T y_intercept() const;
        std::string str() const;
        static Varying<T> growing_varying_with_zero_distance_at_time(T time);
        template<typename U>
        friend std::ostream &operator<<(std::ostream &os, Varying<U> varying);
        bool operator==(Varying<T> rhs) const;
        bool operator!=(Varying<T> rhs) const;
        Varying<T> operator+(T offset) const;
        Varying<T> operator-(T offset) const;
    };

    using Varying32 = Varying<int32_t>;
    using Varying64 = Varying<int64_t>;

    #include "varying.inl"

}

#endif //PYMATCHING2_VARYING_H