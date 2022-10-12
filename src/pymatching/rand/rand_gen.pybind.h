#ifndef PYMATCHING2_RAND_GEN_PYBIND_H
#define PYMATCHING2_RAND_GEN_PYBIND_H

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;
using namespace py::literals;

namespace pm_pybind {

void pybind_rand_gen_methods(py::module &m);

}  // namespace pm_pybind

#endif  // PYMATCHING2_RAND_GEN_PYBIND_H
