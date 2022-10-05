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
        src/pymatching/fill_match/driver/namespaced_main.cc
        src/pymatching/fill_match/driver/io.cc
        src/pymatching/fill_match/driver/mwpm_decoding.cc
        src/pymatching/fill_match/flooder/graph.cc
        src/pymatching/fill_match/flooder/detector_node.cc
        src/pymatching/fill_match/flooder_matcher_interop/compressed_edge.cc
        src/pymatching/fill_match/flooder/graph_fill_region.cc
        src/pymatching/fill_match/flooder/match.cc
        src/pymatching/fill_match/flooder/graph_flooder.cc
        src/pymatching/fill_match/matcher/alternating_tree.cc
        src/pymatching/fill_match/matcher/mwpm.cc
        src/pymatching/fill_match/flooder_matcher_interop/region_edge.cc
        src/pymatching/fill_match/flooder_matcher_interop/mwpm_event.cc
        src/pymatching/fill_match/tracker/flood_check_event.cc
        src/pymatching/fill_match/diagram/animation_main.cc
        src/pymatching/fill_match/diagram/mwpm_diagram.cc
        src/pymatching/fill_match/search/search_graph.cc
        src/pymatching/fill_match/search/search_detector_node.cc
        src/pymatching/fill_match/search/search_flooder.cc
        src/pymatching/fill_match/driver/user_graph.cc
        src/pymatching/rand/rand_gen.cc
        )

set(TEST_FILES
        src/pymatching/fill_match/driver/namespaced_main.test.cc
        src/pymatching/fill_match/driver/io.test.cc
        src/pymatching/fill_match/driver/mwpm_decoding.test.cc
        src/pymatching/fill_match/flooder_matcher_interop/varying.test.cc
        src/pymatching/fill_match/flooder/graph.test.cc
        src/pymatching/fill_match/flooder/detector_node.test.cc
        src/pymatching/fill_match/flooder_matcher_interop/compressed_edge.test.cc
        src/pymatching/fill_match/flooder/graph_fill_region.test.cc
        src/pymatching/fill_match/flooder/match.test.cc
        src/pymatching/fill_match/flooder/graph_flooder.test.cc
        src/pymatching/fill_match/matcher/alternating_tree.test.cc
        src/pymatching/fill_match/matcher/mwpm.test.cc
        src/pymatching/fill_match/tracker/flood_check_event.test.cc
        src/pymatching/fill_match/tracker/radix_heap_queue.test.cc
        src/pymatching/fill_match/flooder_matcher_interop/mwpm_event.test.cc
        src/pymatching/fill_match/tracker/queued_event_tracker.test.cc
        src/pymatching/fill_match/tracker/cyclic.test.cc
        src/pymatching/fill_match/diagram/mwpm_diagram.test.cc
        src/pymatching/fill_match/search/search_graph.test.cc
        src/pymatching/fill_match/search/search_flooder.test.cc
        src/pymatching/fill_match/driver/user_graph.test.cc
        )

set(PERF_FILES
        src/pymatching/perf/main.perf.cc
        src/pymatching/perf/util.perf.cc
        src/pymatching/fill_match/driver/mwpm_decoding.perf.cc
        src/pymatching/fill_match/driver/io.perf.cc
        src/pymatching/fill_match/flooder_matcher_interop/varying.perf.cc
        src/pymatching/fill_match/tracker/radix_heap_queue.perf.cc
        )

set(PYTHON_API_FILES
        src/pymatching/fill_match/driver/user_graph.pybind.cc
        src/pymatching/pymatching.pybind.cc
        )

add_executable(pymatching src/main.cc ${SOURCE_FILES_NO_MAIN})

target_compile_options(pymatching PRIVATE -O3 -march=native)
target_link_options(pymatching PRIVATE -pthread -O3 -march=native)
target_link_libraries(pymatching libstim)

enable_testing()

add_executable(pymatching_tests ${SOURCE_FILES_NO_MAIN} ${TEST_FILES})
target_compile_options(pymatching_tests PRIVATE -fsanitize=address -g -fsanitize=undefined -fno-omit-frame-pointer -coverage)
target_link_options(pymatching_tests PRIVATE -pthread -fsanitize=address -fsanitize=undefined -coverage)
target_link_libraries(pymatching_tests GTest::gtest_main GTest::gmock_main libstim)

add_executable(pymatching_perf ${SOURCE_FILES_NO_MAIN} ${PERF_FILES})
target_compile_options(pymatching_perf PRIVATE -O3 -march=native)
target_link_options(pymatching_perf PRIVATE -pthread -O3)
target_link_libraries(pymatching_perf libstim)

include(GoogleTest)
gtest_discover_tests(pymatching_tests)

if(SKBUILD)
    # Scikit-Build does not add your site-packages to the search path
    # automatically, so we need to add it _or_ the pybind11 specific directory
    # here.
    execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -c
            "import pybind11; print(pybind11.get_cmake_dir())"
            OUTPUT_VARIABLE _tmp_dir
            OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
    list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
endif()

find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(pymatching_python_bindings ${PYTHON_API_FILES} ${SOURCE_FILES_NO_MAIN})
target_link_libraries(pymatching_python_bindings PRIVATE libstim)
set_target_properties(pymatching_python_bindings PROPERTIES OUTPUT_NAME _cpp_pymatching)
add_compile_definitions(PYMATCHING_PYBIND11_MODULE_NAME=_cpp_pymatching)
install(TARGETS pymatching_python_bindings DESTINATION .)