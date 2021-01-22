cmake_minimum_required(VERSION 2.8.12)
project(mwpm)

if(NOT DEFINED CMAKE_MODULE_LINKER_FLAGS_MAINTAINER)
    set(
        CMAKE_MODULE_LINKER_FLAGS_MAINTAINER ""
        CACHE STRING ""
    )
endif()

file(GLOB SOURCES "src/pymatching/*.cpp")

add_subdirectory(lib/pybind11)
pybind11_add_module(_cpp_mwpm ${SOURCES})

add_subdirectory(${CMAKE_SOURCE_DIR}/lib/lemon lib/lemon)

FIND_PACKAGE(Doxygen)

target_link_libraries(_cpp_mwpm PRIVATE lemon)

include_directories(${CMAKE_SOURCE_DIR}/lib/boost_1_72_0)

target_include_directories(_cpp_mwpm PRIVATE
                          "${PROJECT_BINARY_DIR}"
                          "${CMAKE_SOURCE_DIR}/lib/lemon"
                          "${CMAKE_BINARY_DIR}/lib/lemon"
                          )
                    
add_compile_options(-O3)

