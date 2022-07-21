#ifndef PYMATCHING2_VARYING_H
#define PYMATCHING2_VARYING_H

#include <cstdint>
#include <string>

// To do: check that you can only divide weights by at most 2
// Template to take arbitrary unsigned integer type

// Only need two fields. Slope can be int8_t
// y intercept and slope
// Probably one field, 32 bits, last 2 bits be slope, then remainder is y intercept.
// Method to get slope, other to get value at time

// Don't pass by const ref here

namespace pm {
    typedef int32_t time_int;

    class Varying{
    public:
        int32_t data;
        int32_t get_y_intercept();
        uint8_t get_slope();
        int32_t get_distance_at_time(int32_t time);
        Varying then_growing_at_time(int32_t time_of_change);
        Varying then_shrinking_at_time(int32_t time_of_change);
        Varying then_frozen_at_time(int32_t time_of_change);
        std::string str() const;
        Varying();
        explicit Varying(const int32_t& data);
        Varying operator+(const int32_t& offset) const;
        Varying operator-(const int32_t& offset) const;
        friend std::ostream &operator<<(std::ostream &os, const Varying &varying);

        bool operator==(const Varying &rhs) const;

        bool operator!=(const Varying &rhs) const;
    };
    // static method
    Varying growing_varying_with_zero_distance_at_time(int32_t time);
}


#endif //PYMATCHING2_VARYING_H
