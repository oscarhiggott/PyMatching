#ifndef PYMATCHING2_CYCLIC_H
#define PYMATCHING2_CYCLIC_H

#include <sstream>

namespace pm {

template <typename T, class = typename std::enable_if<std::is_unsigned<T>::value>::type>
struct cyclic;

template <typename T>
std::ostream &operator<<(std::ostream &out, const cyclic<T> &c);

/// Helper class for performing comparisons on cyclic spaces.
///
/// A cyclic integer A is LESS THAN a cyclic integer B IFF it takes fewer increments to reach
/// B starting from A than to reach A starting from B.
template <typename T, class>
struct cyclic {
    T value;
    static constexpr T HALF = (T)(((T)-1 >> 1) + (T)1);

    cyclic() : value(0) {
    }

    template <typename O, class = typename std::enable_if<std::is_integral<O>::value>::type>
    explicit cyclic(O value) : value((T)value) {
    }

    inline cyclic<T> operator+(cyclic<T> other) const {
        return cyclic<T>{(T)(value + other.value)};
    }
    inline cyclic<T> operator-(cyclic<T> other) const {
        return cyclic<T>{(T)(value - other.value)};
    }
    inline cyclic<T> operator++() {
        return cyclic<T>(++value);
    }
    inline const cyclic<T> operator++(int) {
        return cyclic<T>(value++);
    }
    inline cyclic<T> operator--() {
        return cyclic<T>(--value);
    }
    inline const cyclic<T> operator--(int) {
        return cyclic<T>(value--);
    }

    template <typename O, class = typename std::enable_if<std::is_integral<O>::value>::type>
    inline cyclic<T> &operator+=(O other) {
        value += other;
        return *this;
    }
    template <typename O, class = typename std::enable_if<std::is_integral<O>::value>::type>
    inline cyclic<T> &operator-=(O other) {
        value -= other;
        return *this;
    }

    inline bool operator<=(cyclic<T> other) const {
        return (T)(other.value - value) < HALF;
    }
    inline bool operator>=(cyclic<T> other) const {
        return (T)(value - other.value) < HALF;
    }
    inline bool operator<(cyclic<T> other) const {
        return (T)(other.value - value - 1) < (T)(HALF - 1);
    }
    inline bool operator>(cyclic<T> other) const {
        return (T)(value - other.value - 1) < HALF - 1;
    }
    inline bool operator==(cyclic<T> other) const {
        return value == other.value;
    }
    inline bool operator!=(cyclic<T> other) const {
        return value != other.value;
    }
    inline bool operator==(T other) const {
        return value == other;
    }
    inline bool operator!=(T other) const {
        return value != other;
    }

    template <typename O, class = typename std::enable_if<std::is_integral<O>::value>::type>
    inline O widen_from_nearby_reference(O reference) const {
        cyclic<T> cyclic_reference((T)reference);
        O result = reference;
        result += (O)(T)(value - cyclic_reference.value);
        if (cyclic_reference > *this) {
            // Move backwards by one full cycle.
            result -= (O)((O)(T)-1 + (O)1);
        }
        return result;
    }

    std::string str() const {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }
};

template <typename T>
std::ostream &operator<<(std::ostream &out, const cyclic<T> &c) {
    out << c.value;
    return out;
}

}  // namespace pm

#endif  // PYMATCHING2_CYCLIC_H
