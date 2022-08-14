cmake_minimum_required(VERSION 3.13)
project(PyMatching2)

include_directories(src)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY out)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY out)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY out)

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

set(SOURCE_FILES_NO_MAIN
        src/pymatching/namespaced_main.cc
        src/pymatching/events.cc
        src/pymatching/graph.cc
        src/pymatching/compressed_edge.cc
        src/pymatching/alternating_tree.cc
        src/pymatching/graph_fill_region.cc
        src/pymatching/region_edge.cc
        src/pymatching/graph_flooder.cc
        src/pymatching/mwpm.cc
        src/pymatching/stim_io.cc
        src/pymatching/mwpm_decoding.cc
        )

set(TEST_FILES
        src/pymatching/varying.test.cc
        src/pymatching/events.test.cc
        src/pymatching/namespaced_main.test.cc
        src/pymatching/graph.test.cc
        src/pymatching/fixed_length_vector.test.cc
        src/pymatching/compressed_edge.test.cc
        src/pymatching/compressed_edge.cc
        src/pymatching/alternating_tree.test.cc
        src/pymatching/graph_fill_region.test.cc
        src/pymatching/graph_flooder.test.cc
        src/pymatching/mwpm.test.cc
        src/pymatching/stim_io.test.cc
        src/pymatching/mwpm_decoding.test.cc
        src/pymatching/bit_bucket_queue.test.cc
        src/pymatching/cyclic.test.cc
        )

set(PERF_FILES
        src/pymatching/perf/main.perf.cc
        src/pymatching/mwpm_decoding.perf.cc
        src/pymatching/perf/util.perf.cc
        src/pymatching/varying.perf.cc
        src/pymatching/bit_bucket_queue.perf.cc
        )

add_executable(pymatching src/main.cc ${SOURCE_FILES_NO_MAIN})

target_compile_options(pymatching PRIVATE -O3)
target_link_options(pymatching PRIVATE -pthread -O3)
target_link_libraries(pymatching libstim)

enable_testing()

add_executable(pymatching_tests ${SOURCE_FILES_NO_MAIN} ${TEST_FILES})
target_compile_options(pymatching_tests PRIVATE -fsanitize=address -g -fsanitize=undefined -fno-omit-frame-pointer)
target_link_options(pymatching_tests PRIVATE -pthread -fsanitize=address -fsanitize=undefined)
target_link_libraries(pymatching_tests GTest::gtest_main libstim)

add_executable(pymatching_perf ${SOURCE_FILES_NO_MAIN} ${PERF_FILES})
target_compile_options(pymatching_perf PRIVATE -O3 -march=native)
target_link_options(pymatching_perf PRIVATE -pthread -O3)
target_link_libraries(pymatching_perf libstim)

include(GoogleTest)
gtest_discover_tests(pymatching_tests)
