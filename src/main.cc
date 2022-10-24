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

#include <iostream>
#include <stdexcept>

#include "pymatching/fill_match.h"

int main(int argc, const char** argv) {
    try {
        return pm::main(argc, argv);
    } catch (std::invalid_argument& ex) {
        std::cerr << ex.what() << "\n";
        return EXIT_FAILURE;
    }
}
