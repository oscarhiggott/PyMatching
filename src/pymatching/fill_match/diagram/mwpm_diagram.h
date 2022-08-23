#ifndef PYMATCHING2_DIAGRAM_MWPM_DIAGRAM_H
#define PYMATCHING2_DIAGRAM_MWPM_DIAGRAM_H

#include "stim.h"
#include "pymatching/fill_match/matcher/mwpm.h"

namespace pm {

std::vector<std::pair<float, float>> dem_detector_coords(const stim::DetectorErrorModel &dem);
void write_mwpm_svg_diagram(const std::vector<std::pair<float, float>> coords, const pm::Mwpm &mwpm, MwpmEvent focused_event, std::ostream &out);

}  // namespace pm

#endif
