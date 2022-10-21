cmake_minimum_required(VERSION 3.4...3.18)
project(PyMatching2)

include_directories(src)
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(SIMD_WIDTH 128)
if (NOT(MSVC))
    set(ARCH_OPT "-O3" "-mno-avx2")
else ()
    set(ARCH_OPT "-O2")
    # https://stackoverflow.com/a/8591946 (RTC1 incompatible with O2 flag)
    STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif ()

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
if (NOT stim_POPULATED)
    FetchContent_Populate(stim)
    add_subdirectory(${stim_SOURCE_DIR})
endif ()

if (NOT (MSVC))
    target_compile_options(libstim PRIVATE -fno-strict-aliasing -fPIC ${ARCH_OPT})
else ()
    target_compile_options(libstim PRIVATE -fPIC ${ARCH_OPT})
endif ()

set(SOURCE_FILES_NO_MAIN
        src/pymatching/sparse_blossom/driver/namespaced_main.cc
        src/pymatching/sparse_blossom/driver/io.cc
        src/pymatching/sparse_blossom/driver/mwpm_decoding.cc
        src/pymatching/sparse_blossom/flooder/graph.cc
        src/pymatching/sparse_blossom/flooder/detector_node.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.cc
        src/pymatching/sparse_blossom/flooder/graph_fill_region.cc
        src/pymatching/sparse_blossom/flooder/match.cc
        src/pymatching/sparse_blossom/flooder/graph_flooder.cc
        src/pymatching/sparse_blossom/matcher/alternating_tree.cc
        src/pymatching/sparse_blossom/matcher/mwpm.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/region_edge.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.cc
        src/pymatching/sparse_blossom/tracker/flood_check_event.cc
        src/pymatching/sparse_blossom/diagram/animation_main.cc
        src/pymatching/sparse_blossom/diagram/mwpm_diagram.cc
        src/pymatching/sparse_blossom/search/search_graph.cc
        src/pymatching/sparse_blossom/search/search_detector_node.cc
        src/pymatching/sparse_blossom/search/search_flooder.cc
        src/pymatching/sparse_blossom/driver/user_graph.cc
        src/pymatching/rand/rand_gen.cc
        )

set(TEST_FILES
        src/pymatching/sparse_blossom/driver/namespaced_main.test.cc
        src/pymatching/sparse_blossom/driver/io.test.cc
        src/pymatching/sparse_blossom/driver/mwpm_decoding.test.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/varying.test.cc
        src/pymatching/sparse_blossom/flooder/graph.test.cc
        src/pymatching/sparse_blossom/flooder/detector_node.test.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.test.cc
        src/pymatching/sparse_blossom/flooder/graph_fill_region.test.cc
        src/pymatching/sparse_blossom/flooder/match.test.cc
        src/pymatching/sparse_blossom/flooder/graph_flooder.test.cc
        src/pymatching/sparse_blossom/matcher/alternating_tree.test.cc
        src/pymatching/sparse_blossom/matcher/mwpm.test.cc
        src/pymatching/sparse_blossom/tracker/flood_check_event.test.cc
        src/pymatching/sparse_blossom/tracker/radix_heap_queue.test.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.test.cc
        src/pymatching/sparse_blossom/tracker/queued_event_tracker.test.cc
        src/pymatching/sparse_blossom/tracker/cyclic.test.cc
        src/pymatching/sparse_blossom/diagram/mwpm_diagram.test.cc
        src/pymatching/sparse_blossom/search/search_graph.test.cc
        src/pymatching/sparse_blossom/search/search_flooder.test.cc
        src/pymatching/sparse_blossom/driver/user_graph.test.cc
        )

set(PERF_FILES
        src/pymatching/perf/main.perf.cc
        src/pymatching/perf/util.perf.cc
        src/pymatching/sparse_blossom/driver/mwpm_decoding.perf.cc
        src/pymatching/sparse_blossom/driver/io.perf.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/varying.perf.cc
        src/pymatching/sparse_blossom/tracker/radix_heap_queue.perf.cc
        )

set(PYTHON_API_FILES
        src/pymatching/sparse_blossom/driver/user_graph.pybind.cc
        src/pymatching/rand/rand_gen.pybind.cc
        src/pymatching/pymatching.pybind.cc
        )

add_executable(pymatching src/main.cc ${SOURCE_FILES_NO_MAIN})

target_compile_options(pymatching PRIVATE ${ARCH_OPT})
if (NOT (MSVC))
    target_link_options(pymatching PRIVATE -pthread ${ARCH_OPT})
endif ()
target_link_libraries(pymatching libstim)

enable_testing()

add_executable(pymatching_tests ${SOURCE_FILES_NO_MAIN} ${TEST_FILES})
target_compile_options(pymatching_tests PRIVATE -fsanitize=address -g -fsanitize=undefined -fno-omit-frame-pointer -coverage)
if (NOT (MSVC))
    target_link_options(pymatching_tests PRIVATE -pthread -fsanitize=address -fsanitize=undefined -coverage)
else ()
    target_link_options(pymatching_tests PRIVATE -fsanitize=address -fsanitize=undefined -coverage)
endif ()
target_link_libraries(pymatching_tests GTest::gtest_main GTest::gmock_main libstim)

add_executable(pymatching_perf ${SOURCE_FILES_NO_MAIN} ${PERF_FILES})
target_compile_options(pymatching_perf PRIVATE ${ARCH_OPT})
if (NOT (MSVC))
    target_link_options(pymatching_perf PRIVATE -pthread -O3)
endif ()
target_link_libraries(pymatching_perf libstim)

include(GoogleTest)
gtest_discover_tests(pymatching_tests)
add_subdirectory(pybind11)
pybind11_add_module(_cpp_pymatching ${PYTHON_API_FILES} ${SOURCE_FILES_NO_MAIN})
target_link_libraries(_cpp_pymatching PRIVATE libstim)
target_compile_options(_cpp_pymatching PRIVATE ${ARCH_OPT})
# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
#target_compile_definitions(_cpp_pymatching
#        PRIVATE VERSION_INFO=${PYMATCHING_VERSION_INFO})