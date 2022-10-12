#ifndef PYMATCHING2_DIAGRAM_MAIN_H
#define PYMATCHING2_DIAGRAM_MAIN_H

namespace pm {

/// Renders image frames as SVGs as matching progresses.
///
/// These frames can then be turned into a video by a tool such as ffmpeg.
/// Here is an example of how to create an animation, assuming stim and pymatching are installed:
///
///     stim gen --rounds=10 --distance=30 --after_clifford_depolarization=0.1 --code repetition_code --task memory >
///     test.stim stim analyze_errors --decompose_errors --fold_loops --in test.stim > test.dem stim detect --shots=1
///     --in test.stim --out_format r8 > test.r8
///
///     pymatching animate \
///         --dem_in test.dem \
///         --dets_in test.r8 \
///         --dets_in_format r8 \
///         --held_frames_per_event 1 \
///         --held_frames_at_start 10 \
///         --held_frames_at_end 10 \
///         --max_growth_between_frames 25 \
///         --max_edge_weight 1000 \
///         --pixels_per_unit_length 20 \
///         --out_dir out_frames
///
///     ffmpeg \
///         -framerate 40 \
///         -pattern_type glob \
///         -i 'out_frames/*.svg' \
///         -vf scale=512:-1 \
///         -c:v libx264 \
///         -vf format=yuv420p \
///         video.mp4
int main_animation(int argc, const char **argv);

}  // namespace pm

#endif
