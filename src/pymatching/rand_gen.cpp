// Copyright 2020 Oscar Higgott

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//      http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <random>
#include "rand_gen.h"

std::mt19937 & global_urng( ) {
static std::mt19937 u{}; return u;
}

void randomize( ) {
static std::random_device rd{};
global_urng().seed( rd() );
}

void set_seed( unsigned seed) {
    global_urng().seed(seed);
}

double rand_float(double from, double to) {
    static std::uniform_real_distribution<> d{};
    using parm_t = decltype(d)::param_type;
    return d( global_urng(), parm_t{from, to} ); 
}