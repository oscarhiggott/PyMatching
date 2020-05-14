#pragma once

#include <random>

std::mt19937 & global_urng();

void randomize();

void set_seed(unsigned s);

double rand_float(double from, double to);