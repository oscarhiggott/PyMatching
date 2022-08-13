#ifndef PYMATCHING2_FIXED_LENGTH_VECTOR_H
#define PYMATCHING2_FIXED_LENGTH_VECTOR_H

#include <cstddef>

namespace pm {

template <class T>
class Vec {
    virtual void push_back(const T& element) = 0;
    virtual std::size_t size() = 0;
};

template <class T, std::size_t N>
class FixedLengthVector : Vec<T> {
    T arr[N];
    size_t num_elements;

   public:
    FixedLengthVector() : num_elements(0) {
    }
    void push_back(const T& element);
    std::size_t size();
    T& operator[](std::size_t index);
};

template <class T, std::size_t N>
inline void FixedLengthVector<T, N>::push_back(const T& element) {
    arr[num_elements] = element;
    num_elements += 1;
};

template <class T, std::size_t N>
inline std::size_t FixedLengthVector<T, N>::size() {
    return num_elements;
};

template <class T, std::size_t N>
inline T& FixedLengthVector<T, N>::operator[](std::size_t index) {
    return arr[index];
}

template <class T>
using Vec16 = FixedLengthVector<T, 16>;

template <class T>
using Vec24 = FixedLengthVector<T, 24>;

template <class T>
using Vector = std::vector<T>;

}  // namespace pm

#endif  // PYMATCHING2_FIXED_LENGTH_VECTOR_H
