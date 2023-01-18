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


Then, using stim to generate some samples, the time per shot 
in microseconds for pymatching 2 was measured (on an M1 Max processor) by running `pymatching.Matching.decode_batch` on 
at least 10000 shots. For example, the number of microseconds per shot can be measured using the following function:
```python
import time
import stim
import pymatching


def time_surface_code_circuit(distance: int, p: float, num_shots: int = 10000) -> float:
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_x",
        distance=distance,
        rounds=distance,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p
    )
    dem = circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    sampler = circuit.compile_detector_sampler()
    shots, actual_observables = sampler.sample(shots=num_shots, separate_observables=True)
    # Decode one shot first to ensure internal C++ representation of the matching graph is fully cached
    matching.decode_batch(shots[0:1, :])
    # Now time decoding the batch
    t0 = time.time()
    matching.decode_batch(shots)
    t1 = time.time()
    microseconds_per_shot = 1e6*(t1-t0)/num_shots
    return microseconds_per_shot
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

If you've looked at the benchmarks in this repository before, you may have noticed that there have been two previous 
versions of the data. In the first version, 
only the X basis was decoded, which did not fully represent the work required to decode a surface code at scale, 
as described above. In the second version, both bases were decoded (and since the problems became 2x bigger, 
the time per round also doubled). However, for both the first and second versions, the timing data was collected by 
decoding shot data from file using the pymatching command line interface. At low p (e.g. around 0.1%), it turned out 
that almost half the time was spent reading the shot data from file. So in the current version the shot data is 
decoded in a batch from memory (see above), with both bases still decoded.
