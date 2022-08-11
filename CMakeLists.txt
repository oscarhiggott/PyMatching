cmake_minimum_required(VERSION 3.15.3)
project(PyMatching2)

set(CMAKE_CXX_STANDARD 20)

#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

FetchContent_MakeAvailable(googletest)

FetchContent_Declare(stim
        GIT_REPOSITORY https://github.com/quantumlib/stim.git
        GIT_TAG v1.9.0)
FetchContent_GetProperties(stim)
if(NOT stim_POPULATED)
    FetchContent_Populate(stim)
    add_subdirectory(${stim_SOURCE_DIR})
endif()

add_executable(
        PyMatching2
        main.cc
        namespaced_main.cc
        events.cc
        graph.cc
        compressed_edge.cc
        alternating_tree.cc
        graph_fill_region.cc
        region_edge.cc
        graph_flooder.cc
        mwpm.cc
        stim_io.cc
        mwpm_decoding.cc
        )

target_compile_options(PyMatching2 PRIVATE -O3)
target_link_options(
        PyMatching2
        PRIVATE
        -pthread
        -O3
)
target_link_libraries(
        PyMatching2
        libstim
)

enable_testing()

add_executable(
        pymatching_tests
        varying.test.cc
        events.cc
        events.test.cc
        namespaced_main.cc
        namespaced_main.test.cc
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
        mwpm.cc
        mwpm.test.cc
        stim_io.cc
        stim_io.test.cc
        mwpm_decoding.cc
)
target_link_options(
        pymatching_tests
        PRIVATE
        -pthread
)
target_link_libraries(
        pymatching_tests
        GTest::gtest_main
        libstim
)

include(GoogleTest)
gtest_discover_tests(pymatching_tests)
