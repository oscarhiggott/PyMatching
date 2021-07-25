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

#include "weighted_stabiliser_graph.h"
#include "stabiliser_graph.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <exception>


struct BlossomFailureException : public std::exception {
    const char * what() const throw();
};


/**
 * @brief A struct containing the output of the minimum weight perfect matching decoder.
 * Contains the correction corresponding to the solution, as well as the total weight 
 * of the solution (with the latter set to -1 if not requested).
 * 
 */
struct MatchingResult {
    /**
     * @brief The correction operator corresponding to the minimum-weight perfect matching.
     * correction[i] is 1 if the ith qubit is flipped and correction[i] is 0 otherwise.
     * 
     */
    py::array_t<std::uint8_t> correction;
    /**
     * @brief The total weight of the edges in the minimum-weight perfect matching.
     * If the weight is not requested by the decoder (return_weight=false), then 
     * weight=-1.
     * 
     */
    double weight;
};

/**
 * @brief Given a stabiliser graph sg and a vector `defects` of indices of nodes that have a -1 syndrome, 
 * find the find the minimum weight perfect matching in the complete graph with nodes in the defects 
 * list, and where the edge between node i and j is given by the distance between i and j in sg. The 
 * distances and shortest paths between nodes in the stabiliser graph sg are all precomputed and this 
 * method returns the exact minimum-weight perfect matching. As a result it is suitable for matching graphs 
 * with a few thousand nodes or less, but will be very memory and compute intensive for larger matching graphs.
 * Returns a noise vector N for which N[i]=1 if qubit_id appeared an odd number of times in the minimum weight 
 * perfect matching and N[i]=0 otherwise.
 * 
 * @param sg A stabiliser graph
 * @param defects The indices of nodes that are associated with a -1 syndrome
 * @return MatchingResult A struct containing the correction vector for the minimum-weight perfect matching and the matching weight. 
 * The matching weight is set to -1 if it is not requested.
 */
MatchingResult LemonDecode(WeightedMatchingGraph& sg, const py::array_t<int>& defects, bool return_weight=false);
/**
 * @brief Given a stabiliser graph `sg`, a vector `defects` of indices of nodes that have a -1 syndrome and 
 * a chosen `num_neighbours`, find the minimum weight perfect matching in the graph V where each defect node 
 * is connected by an edge to each of the `num_neighbours` nearest other defect nodes in sg, and where the 
 * weight of each edge is the distance between the two defect nodes in `sg`.
 * Returns a noise vector N for which N[i]=1 if qubit_id appeared an odd number of times in the minimum weight 
 * perfect matching and N[i]=0 otherwise.
 * 
 * @param sg A stabiliser graph
 * @param defects The indices of nodes that are associated with a -1 syndrome
 * @param num_neighbours The number of closest defects to connect each defect to in the matching graph
 * @return MatchingResult A struct containing the correction vector for the minimum-weight perfect matching and the matching weight. 
 * The matching weight is set to -1 if it is not requested.
 */
MatchingResult LemonDecodeMatchNeighbourhood(WeightedMatchingGraph& sg, std::set<int>& defects,
                                                        int num_neighbours=30, bool return_weight=false);

MatchingResult LocalMatching(
    WeightedMatchingGraph& sg,
    const py::array_t<int>& defects,
    int num_neighbours=30,
    bool return_weight=false,
    int max_attempts=10
    );