# Building command line tool with cmake

```bash
cmake .
make pymatching
# result is at out/pymatching
```

# Predicting observable flips with command line tool

Generate a problem with stim:

```bash
stim gen \
    --rounds=5 \
    --distance=5 \
    --after_clifford_depolarization=0.003 \
    --code surface_code \
    --task rotated_memory_x \
    > circuit.stim
stim analyze_errors \
    --decompose_errors \
    --fold_loops \
    --in circuit.stim \
    > error_model.dem
stim detect \
    --in circuit.stim \
    --shots 10000 \
    --obs_out actual_obs_flips.01 \
    --obs_out_format 01 \
    --out detection_events.b8 \
    --out_format b8
```

Solve problem with pymatching:

```bash
pymatching predict \
    --dem error_model.dem \
    --in detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips.01 \
    --out_format 01
```

Check how well it did:

```bash
echo correct predictions:
paste -d " " predicted_obs_flips.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l
```

# Creating animations with command line tool

```bash
pymatching animate \
    --dem_in test.dem \
    --dets_in test.r8 \
    --dets_in_format r8 \
    --held_frames_per_event 1 \
    --held_frames_at_start 10 \
    --held_frames_at_end 10 \
    --max_growth_between_frames 25 \
    --max_edge_weight 1000 \
    --pixels_per_unit_length 20 \
    --out_dir out_frames
```

# Running C++ tests with cmake

```bash
cmake .
make pymatching_tests
out/pymatching_tests
```

# Running performance profiling with cmake

```bash
cmake .
make pymatching_perf
out/pymatching_perf
```

Use `--target_seconds=#` to change how long each benchmark runs for.

Use `--only=name` to only run one benchmark.

Use `--only=prefix*` to only run benchmarks beginning with a prefix.

# Building the python package locally

```bash
pip install setuptools>=42 pybind11~=2.9.2 cmake>=3.22 scikit-build>=0.15.0
python setup.py bdist
```
