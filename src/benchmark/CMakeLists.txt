add_library(libbenchmark
    util.cc
    mwpm_decoding.cc
    io.cc
    varying.cc
    radix_heap_queue.cc
)

target_link_libraries(libbenchmark libblossom)
