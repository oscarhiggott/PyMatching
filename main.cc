#include <iostream>
#include "varying.perf.h"


int main() {
    auto out = pm::time_get_distance(pm::DIRECT);
    std::cout << std::get<1>(out) << std::endl;
    std::cout << std::get<0>(out) << std::endl;
}
