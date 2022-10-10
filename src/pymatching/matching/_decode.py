from typing import Union, List, Tuple, TYPE_CHECKING
if TYPE_CHECKING:
    import pymatching

import numpy as np


def decode(self: 'pymatching.Matching',
           z: Union[np.ndarray, List[int]],
           *,
           return_weight: bool = False,
           **kwargs
           ) -> Union[np.ndarray, Tuple[np.ndarray, int]]:
    """Decode the syndrome `z` using minimum-weight perfect matching

    Parameters
    ----------
    z : numpy.ndarray
        A binary syndrome vector to decode. The number of elements in
        `z` should equal the number of nodes in the matching graph. If
        `z` is a 1D array, then `z[i]` is the syndrome at node `i` of
        the matching graph. If `z` is 2D then `z[i,j]` is the difference
        (modulo 2) between the (noisy) measurement of stabiliser `i` in time
        step `j+1` and time step `j` (for the case where the matching graph is
        constructed from a check matrix with `repetitions>1`).
    return_weight : bool, optional
        If `return_weight==True`, the sum of the weights of the edges in the
        minimum weight perfect matching is also returned. By default False
    Returns
    -------
    correction : numpy.ndarray or list[int]
        A 1D numpy array of ints giving the minimum-weight correction operator as a
        binary vector. The number of elements in `correction` is one greater than
        the largest fault ID. The ith element of `correction` is 1 if the
        minimum-weight perfect matching (MWPM) found by PyMatching contains an odd
        number of edges that have `i` as one of the `fault_ids`, and is 0 otherwise.
        If each edge in the matching graph is assigned a unique integer in its
        `fault_ids` attribute, then the locations of nonzero entries in `correction`
        correspond to the edges in the MWPM. However, `fault_ids` can instead be used,
        for example, to store IDs of the physical or logical frame changes that occur
        when an edge flips (see the documentation for ``Matching.add_edge`` for more information).
    weight : float
        Present only if `return_weight==True`.
        The sum of the weights of the edges in the minimum-weight perfect
        matching.
    Raises
    ------
    ValueError
        If there is no error consistent with the provided syndrome. Occurs if the syndrome has odd parity in the
        support of a connected component without a boundary.
    Examples
    --------
    >>> import pymatching
    >>> import numpy as np
    >>> H = np.array([[1, 1, 0, 0, 0],
    ...               [0, 1, 1, 0, 0],
    ...               [0, 0, 1, 1, 0],
    ...               [0, 0, 0, 1, 1]])
    >>> m = pymatching.Matching(H)
    >>> z = np.array([0, 1, 0, 0])
    >>> m.decode(z)
    array([1, 1, 0, 0, 0], dtype=uint8)

    Each bit in the correction provided by Matching.decode corresponds to a
    fault_ids. The index of a bit in a correction corresponds to its fault_ids.
    For example, here an error on edge (0, 1) flips fault_ids 2 and 3, as
    inferred by the minimum-weight correction:
    >>> import pymatching
    >>> m = pymatching.Matching()
    >>> m.add_edge(0, 1, fault_ids={2, 3})
    >>> m.add_edge(1, 2, fault_ids=1)
    >>> m.add_edge(2, 0, fault_ids=0)
    >>> m.decode([1, 1, 0])
    array([0, 0, 1, 1], dtype=uint8)

    To decode with a phenomenological noise model (qubits and measurements both suffering
    bit-flip errors), you can provide a check matrix and number of syndrome repetitions to
    construct a matching graph with a time dimension (where nodes in consecutive time steps
    are connected by an edge), and then decode with a 2D syndrome
    (dimension 0 is space, dimension 1 is time):
    >>> import pymatching
    >>> import numpy as np
    >>> np.random.seed(0)
    >>> H = np.array([[1, 1, 0, 0],
    ...               [0, 1, 1, 0],
    ...               [0, 0, 1, 1]])
    >>> m = pymatching.Matching(H, repetitions=5)
    >>> data_qubit_noise = (np.random.rand(4, 5) < 0.1).astype(np.uint8)
    >>> print(data_qubit_noise)
    [[0 0 0 0 0]
     [0 0 0 0 0]
     [0 0 0 0 1]
     [1 1 0 0 0]]
    >>> cumulative_noise = (np.cumsum(data_qubit_noise, 1) % 2).astype(np.uint8)
    >>> syndrome = H@cumulative_noise % 2
    >>> print(syndrome)
    [[0 0 0 0 0]
     [0 0 0 0 1]
     [1 0 0 0 1]]
    >>> syndrome[:,:-1] ^= (np.random.rand(3, 4) < 0.1).astype(np.uint8)
    >>> # Take the parity of consecutive timesteps to construct a difference syndrome:
    >>> syndrome[:,1:] = syndrome[:,:-1] ^ syndrome[:,1:]
    >>> m.decode(syndrome)
    array([0, 0, 1, 0], dtype=uint8)
    """
    try:
        z = np.array(z, dtype=np.uint8)
    except:
        raise TypeError("Syndrome must be of type numpy.ndarray or " \
                        "convertible to numpy.ndarray, not {}".format(z))
    if len(z.shape) == 1 and (self.num_detectors <= z.shape[0]
                              <= self.num_detectors + len(self.boundary)):
        defects = z.nonzero()[0]
    elif len(z.shape) == 2 and z.shape[0] * z.shape[1] == self.num_detectors:
        times, checks = z.T.nonzero()
        defects = times * z.shape[0] + checks
    else:
        raise ValueError("The shape ({}) of the syndrome vector z is not valid.".format(z.shape))
    correction, weight = self._matching_graph.decode(defects)
    if return_weight:
        return correction, weight
    else:
        return correction
