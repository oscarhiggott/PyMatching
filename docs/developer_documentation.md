# Index

- [Build the pymatching command line tool with cmake](#build-cmake)
- [Decode a problem using the pymatching command line tool](#decode)
- [Create a decoding animation using the pymatching command line tool](#decode-animate)
- [Run C++ unit tests using cmake](#cmake-test)
- [Run C++ performance benchmarks tests using cmake](#cmake-perf)
- [Build and install development version of pymatching python package](#pip-install)

# <a name="build-cmake"></a>Build the pymatching command line tool with cmake

```bash
cmake .
make pymatching
# result is at ./pymatching
```

# <a name="decode"></a>Decode a problem using the pymatching command line tool

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

# <a name="decode-animate"></a>Create a decoding animation using the pymatching command line tool

Generate a problem using stim:

```bash
stim gen \
    --rounds=25 \
    --distance=25 \
    --after_clifford_depolarization=0.05 \
    --code repetition_code \
    --task memory \
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

Produce animation frames of solving the problem using pymatching:

```bash
rm out_frames -rf
mkdir out_frames

pymatching animate \
    --dem_in error_model.dem \
    --dets_in detection_events.b8 \
    --dets_in_format b8 \
    --held_frames_per_event 1 \
    --held_frames_at_start 10 \
    --held_frames_at_end 10 \
    --max_growth_between_frames 25 \
    --pixels_per_unit_length 20 \
    --out_dir out_frames
```

Render the produced frames into an mp4 video using ffmpeg:

```bash
ffmpeg \
    -framerate 40 \
    -pattern_type glob \
    -i 'out_frames/*.svg' \
    -vf scale=1024:-1 \
    -c:v libx264 \
    -vf format=yuv420p \
    output_video.mp4
```

# <a name="cmake-test"></a>Run C++ unit tests using cmake

```bash
cmake .
make pymatching_tests
./pymatching_tests
```

# <a name="cmake-perf"></a>Run C++ performance benchmarks tests using cmake

```bash
cmake .
make pymatching_perf
./pymatching_perf
```

Use `--target_seconds=#` to change how long each benchmark runs for.

Use `--only=name` to only run one benchmark.

Use `--only=prefix*` to only run benchmarks beginning with a prefix.

# <a name="pip-install"></a>Build and install development version of pymatching python package

```bash
pip install -e .
```

# <a name="sphinx"></a>Build the Sphinx documentation

First install the sphinx requirements:

```bash
pip install -r docs/sphinx_docs/requirements.txt
```

You will also need to install the latest version of pymatching, and you may also need to [install pandoc](https://pandoc.org/installing.html).

Then, to build the html sphinx docs, run:
```bash
cd docs/sphinx_docs
make html
```

and view `docs/sphinx_docs/build/html/index.html` in a browser.

