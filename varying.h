#ifndef PYMATCHING2_VARYING_H
#define PYMATCHING2_VARYING_H

#include <cstdint>
#include <string>
#include <sstream>


namespace pm {
    typedef int32_t time_int;

    template<typename T>
    class Varying {
    public:
        T data{};
        Varying();
        explicit Varying(T data);
        T get_distance_at_time(T time) const;
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

    template<typename T>
    inline Varying<T>::Varying() = default;

    template<typename T>
    inline Varying<T>::Varying(T data) : data(data) {}

    template<typename T>
    inline T Varying<T>::get_distance_at_time(T time) const {
        if (data & 1) {
            return (data >> 2) + time;
        } else if (data & 2) {
            return (data >> 2) - time;
        } else {
            return data >> 2;
        }
    }

    template<typename T>
    inline Varying<T> Varying<T>::then_growing_at_time(T time_of_change) const {
        T d = get_distance_at_time(time_of_change);
        return pm::Varying<T>(((d - time_of_change) << 2) + 1);
    }

    template<typename T>
    inline Varying<T> Varying<T>::then_shrinking_at_time(T time_of_change) const {
        T d = get_distance_at_time(time_of_change);
        return pm::Varying<T>(((d + time_of_change) << 2) + 2);
    }

    template<typename T>
    inline Varying<T> Varying<T>::then_frozen_at_time(T time_of_change) const {
        T d = get_distance_at_time(time_of_change);
        return pm::Varying<T>(d << 2);
    }

    template<typename T>
    inline std::string Varying<T>::str() const {
        std::stringstream s;
        s << *this;
        return s.str();
    }

    template<typename T>
    inline Varying<T> Varying<T>::growing_varying_with_zero_distance_at_time(T time) {
        return Varying<T>(((-time) << 2) + 1);
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &os, Varying<T> varying) {
        os << (varying.data >> 2);
        if (varying.data & 1){
            os << " + t";
        } else if (varying.data & 2){
            os << " - t";
        }
        return os;
    }

    template<typename T>
    inline bool Varying<T>::operator==(Varying<T> rhs) const {
        return data == rhs.data;
    }

    template<typename T>
    inline bool Varying<T>::operator!=(Varying<T> rhs) const {
        return !(rhs == *this);
    }

    template<typename T>
    inline Varying<T> Varying<T>::operator+(T offset) const {
        return Varying<T>(this->data + (offset<<2));
    }

    template<typename T>
    inline Varying<T> Varying<T>::operator-(T offset) const {
        return pm::Varying<T>(this->data - (offset << 2));
    }
}

#endif //PYMATCHING2_VARYING_H
