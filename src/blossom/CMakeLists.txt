add_library(libblossom
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
target_link_libraries(libblossom libstim)
