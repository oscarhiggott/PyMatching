#include "varying.h"
#include <sstream>

uint8_t pm::Varying::get_slope() {
    return data & 3;
}

int32_t pm::Varying::get_y_intercept() {
    return data >> 2;
}

int32_t pm::Varying::get_distance_at_time(int32_t time) {
    if (data & 1) {
        return (data >> 2) + time;
    } else if (data & 2) {
        return (data >> 2) - time;
    } else {
        return data >> 2;
    }
}

pm::Varying::Varying() {}

pm::Varying::Varying(const int32_t& data) : data(data) {}

pm::Varying pm::growing_varying_with_zero_distance_at_time(int32_t time) {
    return Varying(((-time) << 2) + 1);
}

pm::Varying pm::Varying::operator+(const int32_t &offset) const {
    return Varying(this->data + (offset<<2));
}

pm::Varying pm::Varying::operator-(const int32_t &offset) const {
    return pm::Varying(this->data - (offset << 2));
}


std::ostream &pm::operator<<(std::ostream &os, const pm::Varying &varying) {
    os << (varying.data >> 2);
    if (varying.data & 1){
        os << " + t";
    } else if (varying.data & 2){
        os << " - t";
    }
    return os;
}

pm::Varying pm::Varying::then_growing_at_time(int32_t time_of_change) {
    int32_t d = get_distance_at_time(time_of_change);
    return pm::Varying(((d - time_of_change) << 2) + 1);
}

pm::Varying pm::Varying::then_shrinking_at_time(int32_t time_of_change) {
    int32_t d = get_distance_at_time(time_of_change);
    return pm::Varying(((d + time_of_change) << 2) + 2);
}

pm::Varying pm::Varying::then_frozen_at_time(int32_t time_of_change) {
    int32_t d = get_distance_at_time(time_of_change);
    return pm::Varying(d << 2);
}

bool pm::Varying::operator==(const pm::Varying &rhs) const {
    return data == rhs.data;
}

bool pm::Varying::operator!=(const pm::Varying &rhs) const {
    return !(rhs == *this);
}

std::string pm::Varying::str() const {
    std::stringstream s;
    s << *this;
    return s.str();
}

