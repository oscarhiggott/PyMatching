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
