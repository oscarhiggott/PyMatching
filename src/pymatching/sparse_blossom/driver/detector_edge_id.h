#ifndef PYMATCHING2_DETECTOR_EDGE_ID_H_
#define PYMATCHING2_DETECTOR_EDGE_ID_H_
#include <cstdint>
#include <string>

namespace stim_util {

/// Identifies an edge or boundary edge in a detector graph.
///
/// The edge is identified by the one or two detectors it flips.
/// Boundary edges store the same detector twice.
/// The detectors are kept sorted to avoid equality ambiguity.
struct DetectorEdgeId {
    /// The detector index (to find detector in appropriate vector) for the first detector.
    uint64_t d1;
    /// The detector index (to find detector in appropriate vector) for the second detector.
    /// Will be same as `d1` if this `DetectorEdgeId` is for a boundary edge.
    uint64_t d2;

    DetectorEdgeId();
    DetectorEdgeId(uint64_t detector_id_1);
    DetectorEdgeId(uint64_t detector_id_1, uint64_t detector_id_2);

    const uint64_t *begin() const;
    const uint64_t *end() const;
    bool is_boundary_edge() const;

    bool operator==(const DetectorEdgeId &other) const;
    bool operator!=(const DetectorEdgeId &other) const;
    bool operator<(const DetectorEdgeId &other) const;
    std::string str() const;
};

}  // namespace stim_util
#endif
