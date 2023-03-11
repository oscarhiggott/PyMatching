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

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

FetchContent_MakeAvailable(googletest)

set(TEST_FILES
    src/blossom/driver/namespaced_main.test.cc
    src/blossom/driver/io.test.cc
    src/blossom/driver/mwpm_decoding.test.cc
    src/blossom/flooder_matcher_interop/varying.test.cc
    src/blossom/flooder/graph.test.cc
    src/blossom/flooder/detector_node.test.cc
    src/blossom/flooder_matcher_interop/compressed_edge.test.cc
    src/blossom/flooder/graph_fill_region.test.cc
    src/blossom/flooder/match.test.cc
    src/blossom/flooder/graph_flooder.test.cc
    src/blossom/matcher/alternating_tree.test.cc
    src/blossom/matcher/mwpm.test.cc
    src/blossom/tracker/flood_check_event.test.cc
    src/blossom/tracker/radix_heap_queue.test.cc
    src/blossom/flooder_matcher_interop/mwpm_event.test.cc
    src/blossom/tracker/queued_event_tracker.test.cc
    src/blossom/tracker/cyclic.test.cc
    src/blossom/diagram/mwpm_diagram.test.cc
    src/blossom/search/search_graph.test.cc
    src/blossom/search/search_flooder.test.cc
    src/blossom/driver/user_graph.test.cc
)

set(PYTHON_API_FILES
    src/blossom/driver/user_graph.pybind.cc
    src/pymatching/rand/rand_gen.pybind.cc
    src/pymatching/pymatching.pybind.cc
)

enable_testing()
add_executable(pymatching_tests ${TEST_FILES})
target_link_libraries(pymatching_tests libblossom)
target_compile_options(pymatching_tests PRIVATE -fsanitize=address -g -fsanitize=undefined -fno-omit-frame-pointer -coverage)
if (NOT (MSVC))
    target_link_options(pymatching_tests PRIVATE -pthread -fsanitize=address -fsanitize=undefined -coverage)
else ()
    target_link_options(pymatching_tests PRIVATE -fsanitize=address -fsanitize=undefined -coverage)
endif ()
target_link_libraries(pymatching_tests GTest::gtest_main GTest::gmock_main libstim)

include(GoogleTest)
gtest_discover_tests(pymatching_tests)
add_subdirectory(pybind11)
pybind11_add_module(_cpp_pymatching ${PYTHON_API_FILES})
target_link_libraries(_cpp_pymatching PRIVATE libblossom)
target_compile_options(_cpp_pymatching PRIVATE ${ARCH_OPT})
# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
#target_compile_definitions(_cpp_pymatching
#        PRIVATE VERSION_INFO=${PYMATCHING_VERSION_INFO})
