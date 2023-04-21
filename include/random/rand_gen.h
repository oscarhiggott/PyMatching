// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PYMATCHING2_RAND_GEN_H
#define PYMATCHING2_RAND_GEN_H

#include <random>

namespace pm {
std::mt19937& global_urng();

void randomize();

/**
 * @brief Set the seed for the mt19937 random number generator
 *
 * @param s
 */
void set_seed(unsigned seed);

/**
 * @brief A random double chosen uniformly at random between `from` and `to`
 *
 * @param from
 * @param to
 * @return double
 */
double rand_float(double from, double to);
}  // namespace pm

#endif  // PYMATCHING2_RAND_GEN_H
