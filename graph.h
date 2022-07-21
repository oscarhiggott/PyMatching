#ifndef PYMATCHING2_GRAPH_H
#define PYMATCHING2_GRAPH_H


#include <vector>

typedef uint64_t obs_int;
typedef uint8_t weight_int;

class GraphFillRegion;
class TentativeNeighborInteractionEvent;

class DetectorNode{
public:
    obs_int observables_crossed_from_source;
    DetectorNode* reached_from_source;
    int distance_from_source; // DO I NEED THIS???
    GraphFillRegion* region_that_arrived;

    // Eventually, make optionally fixed size arrays with templates.
    std::vector<DetectorNode*> neighbors;
    std::vector<weight_int> neighbor_weights;
    std::vector<obs_int> neighbor_observables;
    std::vector<TentativeNeighborInteractionEvent*> neighbor_schedules;
    // no neighbor back index
};


class Graph {
public:
    std::vector<DetectorNode> nodes;

    void add_edge(int u, int v, weight_int weight, obs_int observables);
};


#endif //PYMATCHING2_GRAPH_H
