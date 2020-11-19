#include "weighted_stabiliser_graph.h"
#include "stabiliser_graph.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


int LemonTest();
double LemonMatchingTest();
py::array_t<std::uint8_t> DecodeLemonMatchNeighbourhood(WeightedStabiliserGraph& sg, const py::array_t<int>& defects, int num_neighbours);