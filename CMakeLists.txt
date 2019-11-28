cmake_minimum_required(VERSION 2.8.12)
project(mwpm)

file(GLOB SOURCES "src/mwpm/*.cpp")

add_subdirectory(lib/blossom5-v2.05.src)

add_subdirectory(lib/pybind11)
pybind11_add_module(_cpp_mwpm ${SOURCES})

target_link_libraries(_cpp_mwpm PRIVATE Blossom5)

target_include_directories(_cpp_mwpm PRIVATE
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/lib/blossom5-v2.05.src"
                          )
                    
add_compile_options(-O3)