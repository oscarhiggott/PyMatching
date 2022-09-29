#include "pymatching/fill_match/driver/python_api_graph.h"

double pm::merge_weights(double a, double b){
    auto sgn = std::copysign(1, a) * std::copysign(1, b);
    auto signed_min = sgn * std::min(std::abs(a), std::abs(b));
    return signed_min + std::log(1 + std::exp(-std::abs(a + b))) - std::log(1 + std::exp(-std::abs(a - b)));
}