cmake_minimum_required(VERSION 3.4...3.18)
project(PyMatching2)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SIMD_WIDTH 128)

include_directories(src)
include_directories(include)

add_subdirectory(src/blossom)
add_subdirectory(src/benchmark)

add_subdirectory(apps/pymatching)
add_subdirectory(apps/benchmarks)

add_subdirectory(tests/cpp)


set(PYTHON_API_FILES
    src/blossom/driver/user_graph.pybind.cc
    src/pymatching/rand/rand_gen.pybind.cc
    src/pymatching/pymatching.pybind.cc
)
add_subdirectory(pybind11)
pybind11_add_module(_cpp_pymatching ${PYTHON_API_FILES})
target_link_libraries(_cpp_pymatching PRIVATE libblossom)
target_compile_options(_cpp_pymatching PRIVATE ${ARCH_OPT})
# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
#target_compile_definitions(_cpp_pymatching
#        PRIVATE VERSION_INFO=${PYMATCHING_VERSION_INFO})
