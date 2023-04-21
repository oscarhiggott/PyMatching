enable_testing()

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

FetchContent_MakeAvailable(googletest)

add_executable(pymatching_tests
    diagram_tests/test_mwpm_diagram.cc
    driver_tests/test_io.cc
    driver_tests/test_mwpm_decoding.cc
    driver_tests/test_namespaced_main.cc
    driver_tests/test_user_graph.cc
    flooder_matcher_interop_tests/test_compressed_edge.cc
    flooder_matcher_interop_tests/test_mwpm_event.cc
    flooder_matcher_interop_tests/test_varying.cc
    flooder_tests/test_detector_node.cc
    flooder_tests/test_graph_fill_region.cc
    flooder_tests/test_graph_flooder.cc
    flooder_tests/test_graph.cc
    flooder_tests/test_match.cc
    matcher_tests/test_alternating_tree.cc
    matcher_tests/test_mwpm.cc
    search_tests/test_search_flooder.cc
    search_tests/test_search_graph.cc
    tracker_tests/test_cyclic.cc
    tracker_tests/test_flood_check_event.cc
    tracker_tests/test_queued_event_tracker.cc
    tracker_tests/test_radix_heap_queue.cc
)

target_compile_options(pymatching_tests PRIVATE -fsanitize=address -g -fsanitize=undefined -fno-omit-frame-pointer -coverage)
if (NOT (MSVC))
    target_link_options(pymatching_tests PRIVATE -pthread -fsanitize=address -fsanitize=undefined -coverage)
else ()
    target_link_options(pymatching_tests PRIVATE -fsanitize=address -fsanitize=undefined -coverage)
endif ()
target_link_libraries(pymatching_tests blossom GTest::gtest_main GTest::gmock_main)

include(GoogleTest)
gtest_discover_tests(pymatching_tests)
