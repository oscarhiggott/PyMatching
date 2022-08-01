cmake_minimum_required(VERSION 3.21)
project(PyMatching2)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

FetchContent_MakeAvailable(googletest)

add_executable(PyMatching2 main.cc graph.cc region_edge.cc alternating_tree.cc)

enable_testing()

add_executable(
        pymatching_tests
        varying.test.cc
        events.cc
        events.test.cc
        graph.test.cc
        graph.cc
        fixed_length_vector.test.cc
        compressed_edge.test.cc
        compressed_edge.cc
        alternating_tree.test.cc
        alternating_tree.cc
        graph_fill_region.cc
        graph_fill_region.test.cc
        region_edge.cc
        graph_flooder.cc
        graph_flooder.test.cc
)
target_link_libraries(
        pymatching_tests
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(pymatching_tests)
