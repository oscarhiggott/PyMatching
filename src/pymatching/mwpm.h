#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "stabiliser_graph.h"
#include "weighted_stabiliser_graph.h"

namespace py = pybind11;

py::array_t<std::uint8_t> Decode(const IStabiliserGraph& sg, const py::array_t<int>& defects);
py::array_t<std::uint8_t> DecodeMatchNeighbourhood(WeightedStabiliserGraph& sg, const py::array_t<int>& defects, int num_neighbours);