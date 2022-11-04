# Copyright 2022 PyMatching Contributors

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pymatching

import doctest
import pytest


def test_matching_docstrings():
    pytest.importorskip("stim")
    doctest.testmod(pymatching.matching, raise_on_error=True)


def test_cpp_pymatching_docstrings():
    doctest.testmod(pymatching._cpp_pymatching, raise_on_error=True)
