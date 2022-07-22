#ifndef PYMATCHING2_VARYING_H
#define PYMATCHING2_VARYING_H

#include <cstdint>
#include <string>
#include <sstream>

// To do: check that you can only divide weights by at most 2
// Template to take arbitrary unsigned integer type


namespace pm {
    typedef int32_t time_int;

    class Varying {
    public:
        int32_t data{};
        Varying();
        explicit Varying(const int32_t &data);
        int32_t get_distance_at_time(int32_t time) const;
        Varying then_growing_at_time(int32_t time_of_change) const;
        Varying then_shrinking_at_time(int32_t time_of_change) const;
        Varying then_frozen_at_time(int32_t time_of_change) const;
        std::string str() const;
        static Varying growing_varying_with_zero_distance_at_time(int32_t time);
        friend std::ostream &operator<<(std::ostream &os, Varying varying);
        bool operator==(Varying rhs) const;
        bool operator!=(Varying rhs) const;
        Varying operator+(int32_t offset) const;
        Varying operator-(int32_t offset) const;
    };

    inline Varying::Varying() = default;

    inline Varying::Varying(const int32_t& data) : data(data) {}

    inline int32_t Varying::get_distance_at_time(int32_t time) const {
        if (data & 1) {
            return (data >> 2) + time;
        } else if (data & 2) {
            return (data >> 2) - time;
        } else {
            return data >> 2;
        }
    }

    inline Varying Varying::then_growing_at_time(int32_t time_of_change) const {
        int32_t d = get_distance_at_time(time_of_change);
        return pm::Varying(((d - time_of_change) << 2) + 1);
    }

    inline Varying Varying::then_shrinking_at_time(int32_t time_of_change) const {
        int32_t d = get_distance_at_time(time_of_change);
        return pm::Varying(((d + time_of_change) << 2) + 2);
    }

    inline Varying Varying::then_frozen_at_time(int32_t time_of_change) const {
        int32_t d = get_distance_at_time(time_of_change);
        return pm::Varying(d << 2);
    }

    inline std::string Varying::str() const {
        std::stringstream s;
        s << *this;
        return s.str();
    }

    inline Varying Varying::growing_varying_with_zero_distance_at_time(int32_t time) {
        return Varying(((-time) << 2) + 1);
    }

    inline std::ostream &operator<<(std::ostream &os, Varying varying) {
        os << (varying.data >> 2);
        if (varying.data & 1){
            os << " + t";
        } else if (varying.data & 2){
            os << " - t";
        }
        return os;
    }

    inline bool Varying::operator==(Varying rhs) const {
        return data == rhs.data;
    }

    inline bool Varying::operator!=(Varying rhs) const {
        return !(rhs == *this);
    }

    inline Varying Varying::operator+(int32_t offset) const {
        return Varying(this->data + (offset<<2));
    }

    inline Varying Varying::operator-(int32_t offset) const {
        return pm::Varying(this->data - (offset << 2));
    }
}

#endif //PYMATCHING2_VARYING_H
