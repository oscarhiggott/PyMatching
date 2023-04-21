find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_cpp_pymatching user_graph.cc rand_gen.cc pymatching.cc)
target_link_libraries(_cpp_pymatching PRIVATE blossom)
target_compile_options(_cpp_pymatching PRIVATE ${ARCH_OPT})
# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
#target_compile_definitions(_cpp_pymatching
#        PRIVATE VERSION_INFO=${PYMATCHING_VERSION_INFO})
