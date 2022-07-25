#ifndef PYMATCHING2_VARYING_H
#define PYMATCHING2_VARYING_H

#include <cstdint>
#include <string>
#include <sstream>
#include <stdexcept>


namespace pm {
    typedef int32_t time_int;

    template<typename T>
    class Varying {
    public:
        T data{};
        Varying();
        explicit Varying(T data);
        T get_distance_at_time(T time) const;
        T get_shrinking_distance_at_time(T time) const;
        T get_growing_distance_at_time(T time) const;
        T get_frozen_distance_at_time(T time) const;
        T time_of_x_intercept();
        T time_of_x_intercept_when_added_to(Varying<T> other);
        Varying<T> then_growing_at_time(T time_of_change) const;
        Varying<T> then_shrinking_at_time(T time_of_change) const;
        Varying<T> then_frozen_at_time(T time_of_change) const;
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
