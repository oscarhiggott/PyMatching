include(FetchContent)

FetchContent_Declare(
    stim
    GIT_REPOSITORY https://github.com/quantumlib/stim.git
    GIT_TAG 15f9dd72210d0a7fdc8602d3d03ebb7ce9271d9d
)
FetchContent_MakeAvailable(stim)

if (NOT (MSVC))
    target_compile_options(libstim PRIVATE -fno-strict-aliasing -fPIC ${ARCH_OPT})
else ()
    target_compile_options(libstim PRIVATE -fPIC ${ARCH_OPT})
endif ()

add_library(blossom
    driver/namespaced_main.cc
    driver/io.cc
    driver/mwpm_decoding.cc
    flooder/graph.cc
    flooder/detector_node.cc
    flooder_matcher_interop/compressed_edge.cc
    flooder/graph_fill_region.cc
    flooder/match.cc
    flooder/graph_flooder.cc
    matcher/alternating_tree.cc
    matcher/mwpm.cc
    flooder_matcher_interop/region_edge.cc
    flooder_matcher_interop/mwpm_event.cc
    tracker/flood_check_event.cc
    diagram/animation_main.cc
    diagram/mwpm_diagram.cc
    search/search_graph.cc
    search/search_detector_node.cc
    search/search_flooder.cc
    driver/user_graph.cc
    # Will remove when random generation is a separate library
    ../../src/pymatching/rand/rand_gen.cc
)
target_link_libraries(blossom libstim)
