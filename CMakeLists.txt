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

add_executable(PyMatching2 main.cc graph.cc region_path.cc alternating_tree.cc)

enable_testing()

add_executable(
        varying_test
        varying.test.cc
)
target_link_libraries(
        varying_test
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(varying_test)
