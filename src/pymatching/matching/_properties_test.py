# Copyright 2022 Oscar Higgott

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from pymatching.matching import Matching


def test_set_boundary_nodes():
    m = Matching([[1, 1, 0], [0, 1, 1]])
    m.set_boundary_nodes({1, 2, 4})
    assert m.boundary == {1, 2, 4}
