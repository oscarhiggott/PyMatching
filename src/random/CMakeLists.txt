add_library(random rand_gen.cc)
target_compile_options(random PRIVATE -fPIC)
