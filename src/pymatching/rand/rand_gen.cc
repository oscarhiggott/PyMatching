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

#include "pymatching/rand/rand_gen.h"

#include <random>

std::mt19937& pm::global_urng() {
    static std::mt19937 u{};
    return u;
}

void pm::randomize() {
    static std::random_device rd{};
    pm::global_urng().seed(rd());
}

void pm::set_seed(unsigned seed) {
    pm::global_urng().seed(seed);
}

double pm::rand_float(double from, double to) {
    static std::uniform_real_distribution<> d{};
    using parm_t = decltype(d)::param_type;
    return d(global_urng(), parm_t{from, to});
}