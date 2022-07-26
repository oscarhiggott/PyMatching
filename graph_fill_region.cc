#include "graph_fill_region.h"
#include "varying.h"

pm::GraphFillRegion::GraphFillRegion()
    : blossom_parent(nullptr), alt_tree_node(nullptr), radius((0 << 2) + 1) {}
