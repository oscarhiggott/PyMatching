cmake_minimum_required(VERSION 3.4...3.18)
project(PyMatching2)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SIMD_WIDTH 128)

include_directories(include)

add_subdirectory(src/blossom)
add_subdirectory(src/random)
add_subdirectory(src/benchmark)

if (IS_DIRECTORY ${pybind11_DIR})
    add_subdirectory(src/python_bindings)
endif()

add_subdirectory(apps/pymatching)
add_subdirectory(apps/benchmarks)

add_subdirectory(tests/cpp)
