The stim circuits used for benchmarks in the subdirectories of this directory are surface code circuits with 
circuit-level noise. Each subdirectory contains circuits at for some fixed error rate `p`. 
The stim circuits can be regenerated using the following function:

```python
import stim


def save_benchmark_circuit(distance: int, p: float) -> None:
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_x",
        distance=distance,
        rounds=distance,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p
    )
    circuit.to_file(f"surface_code_rotated_memory_x_{distance}_{p}.stim")
```

for a given error rate `p` and code distance `distance` in `{5, 7, 9, 13, 17, 23, 29, 39, 50}`.

From a stim circuit `circuit`, the corresponding `stim.DetectorErrorModel` used to configure the decoder can be 
generated using `circuit.detector_error_model(decompose_errors=True)`.

Then, using stim to generate some samples in b8 format (in this case with appended observables), the time per shot 
in microseconds for pymatching 2 was measured (on an M1 Max processor) using the pymatching command line tool, where here 
`$samples_fn` is the filename of the b8 samples file and `$dem_fn` is the filename of the detector error model:

```shell
pymatching count_mistakes --in $samples_fn --in_format b8 --dem $dem_fn --in_includes_appended_observables --time 2>&1 >/dev/null | sed -n -E 's/Decoding time per shot: (.+)us/\1/p' 
```

In the figure in each subdirectory, 
the time per shot is divided by the number of rounds (equal to the distance) to give the time per _round_.
The time per _shot_ in microseconds using PyMatching v2.0 is given in `pymatching_v2.csv`.
For comparison, the plot also includes timing data for the same stim circuits using PyMatching v0.7 
(with `num_neighbours=30`) and using NetworkX.

Note that here the detector error models decoded by pymatching contain detectors associated with both X basis _and_ 
Z basis measurements, despite the fact that only the X basis measurements have any impact on the 
logical error rate of the MWPM decoder when measuring the X logical observable as we do here.
However, we choose to include both the X-type and Z-type detectors, as this more accurately represents the work that 
must be done by the decoder within the wider context of a fault-tolerant quantum computation, 
where the basis of the logical measurement is not always known apriori (both X and Z logical operators must be 
preserved) and the X and Z matching graphs can become connected (such as when implementing a logical S gate in the 
surface code by [braiding twist defects](https://arxiv.org/abs/1609.04673)).
