from ._cpp_mwpm import (PerfectMatching, Options, randomize, set_seed,
                        rand_float)
from .matching import *

randomize() # Set random seed using std::random_device
