#ifndef PYMATCHING2_USER_GRAPH_PYBIND_H
#define PYMATCHING2_USER_GRAPH_PYBIND_H

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pymatching/fill_match/driver/user_graph.h"

namespace py = pybind11;

namespace pm_pybind {

struct CompressedSparseColumnCheckMatrix {
    explicit CompressedSparseColumnCheckMatrix(const py::object&matrix);
    py::array_t<uint8_t> data;
    py::array_t<int64_t> indices;
    py::array_t<int64_t> indptr;
    size_t num_rows;
    size_t num_cols;
};

template <typename T>
py::array_t<T> vec_to_array(std::vector<T> *vec) {
    auto err_capsule = py::capsule(vec, [](void *x) {
        delete reinterpret_cast<std::vector<T> *>(x);
    });
    return py::array_t<T>(vec->size(), vec->data(), err_capsule);
}

py::class_<pm::UserGraph> pybind_user_graph(py::module &m);

void pybind_user_graph_methods(py::module &m, py::class_<pm::UserGraph> &g);

}  // namespace pm_pybind

#endif  // PYMATCHING2_USER_GRAPH_PYBIND_H
