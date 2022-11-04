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

#include "pymatching/rand/rand_gen.pybind.h"

#include "pymatching/rand/rand_gen.h"

void pm_pybind::pybind_rand_gen_methods(py::module &m) {
    m.def("randomize", &pm::randomize, R"pbdoc(
        Choose a random seed using std::random_device

        Examples
        --------
            >>> import pymatching
            >>> pymatching.randomize()
     )pbdoc");
    m.def("set_seed", &pm::set_seed, "seed"_a, R"pbdoc(
        Sets the seed of the random number generator

        Parameters
        ----------
        seed: int
            The seed for the random number generator (must be non-negative)

        Examples
        --------
        >>> import pymatching
        >>> pymatching.set_seed(10)

     )pbdoc");

    m.def("rand_float", &pm::rand_float, "from"_a, "to"_a, R"pbdoc(
        Generate a floating point number chosen uniformly at random
        over the interval between `from` and `to`
        Parameters
        ----------
        from: float
            Smallest float that can be drawn from the distribution
        to: float
            Largest float that can be drawn from the distribution
        Returns
        -------
        float
            The random float
     )pbdoc");
}
