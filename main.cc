#include <iostream>
#include <bitset>
#include "graph.h"
#include "events.h"
#include "region_path.h"
#include "graph_fill_region.h"
#include "graph_flooder.h"
#include "mwpm.h"
#include "varying.h"


int main() {
    pm::Varying32 v = pm::Varying32::growing_varying_with_zero_distance_at_time(-10);
    std::cout << v;
    return 0;
}
