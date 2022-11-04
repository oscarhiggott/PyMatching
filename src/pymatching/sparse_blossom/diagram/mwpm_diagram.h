// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PYMATCHING2_DIAGRAM_MWPM_DIAGRAM_H
#define PYMATCHING2_DIAGRAM_MWPM_DIAGRAM_H

#include "pymatching/sparse_blossom/matcher/mwpm.h"
#include "stim.h"

namespace pm {

/// Draws a series of frames showing a decoding progressing from start to finish.
///
/// The frames can be turned into a movie or a GIF by a tool like ffmpeg. For example:
///
///     ffmpeg \
///         -framerate 30 \
///         -pattern_type glob \
///         -i 'WHERE_YOU_WROTE_THE_FRAMES/*.svg' \
///         output_video.mp4
///
/// Args:
///     mwpm: The initialized decoder, ready to accept new detection events.
///     coords: The 2d location to use for each detector node. Must be the correct length.
///     boundary_coords: For each node, where to end a line to the boundary starting from that node.
///         For nodes not beside the boundary, the value is not used. Must be the correct length.
///     detection_events: The problem to decode. The indices of detector nodes that have a
///         detection event.
///     file_path_prefix: Each frame goes to a separate file. This is the common prefix to use for
///         all the paths (e.g. this could be the directory INCLUDING THE SLASH).
///     print_progress: If true, progress messages are printed to stderr.
///     frames_per_mwpm_event: When events such as region collisions and blossom shatters occur,
///         this is how many frames are spent highlighting that event. Set to 0 to not highlight
///         the events (resulting in a smooth animation without pauses where region growth
///         temporarily halts).
///     hold_frames_at_start: The number of frames to show the initial state, with just the
///         detection events showing.
///     hold_frames_at_end: The number of frames to show the final state, with all regions frozen.
///     max_growth_between_frames: Determines how much regions can grow between frames. Setting this
///         to a low number creates smoother animations, but requires outputting more frames. Set to
///         0 to completely disable growth frames.
void write_animated_decoding_svg_frames(
    pm::Mwpm &mwpm,
    const std::vector<std::pair<float, float>> &coords,
    const std::vector<std::pair<float, float>> &boundary_coords,
    const std::vector<uint64_t> detection_events,
    const std::string &file_path_prefix,
    bool print_progress,
    size_t frames_per_mwpm_event,
    size_t hold_frames_at_start,
    size_t hold_frames_at_end,
    size_t max_growth_between_frames);

/// Draws the current state of the decoder as an SVG image.
///
/// Args:
///     coords: The 2d location to use for each detector node. Must be the correct length.
///     boundary_coords: For each node, where to end a line to the boundary starting from that node.
///         For nodes not beside the boundary, the value is not used. Must be the correct length.
///     mwpm: The decoder in a state to be drawn.
///     focused_event: An event to highlight in the drawing. Set to Event::no_event() to not
///         focus anything.
///     out: Where to write the SVG image.
void write_decoder_state_as_svg(
    const std::vector<std::pair<float, float>> &coords,
    const std::vector<std::pair<float, float>> &boundary_coords,
    const pm::Mwpm &mwpm,
    pm::MwpmEvent focused_event,
    std::ostream &out);

/// Picks reasonable 2d coordinates to use when drawing the state of the decoder.
///
/// Args:
///     dem: The detector error model.
///     pixels_per_unit_length: How much to multiply coordinates by before drawing.
///
/// Returns:
///     first item of pair: The projected coordinates for each detector node.
///     second item of pair: The inferred boundary coordinates for each detector node.
std::pair<std::vector<std::pair<float, float>>, std::vector<std::pair<float, float>>> pick_coords_for_drawing_from_dem(
    const stim::DetectorErrorModel &dem, float pixels_per_unit_length);

int main_animated(int argc, const char **argv);

}  // namespace pm

#endif
