#include "stim_io.h"

struct ModelFlattener {
    size_t det_offset;

    void helper(const stim::DetectorErrorModel& model, uint64_t reps) {
        for (size_t i = 0; i < reps; i++) {
            for (auto instruction : model.instructions) {
                switch (instruction.type) {
                    case stim::DEM_REPEAT_BLOCK:
                        helper(model.blocks[instruction.target_data[1].data],
                               instruction.target_data[0].data);
                        break;
                    case stim::DEM_ERROR:
                        for (size_t j = 0; i < instruction.target_data.size(); i++) {
                            if (instruction.target_data[j].is_relative_detector_id()){

                            } else if (instruction.target_data[j].is_observable_id()) {

                            } else if (instruction.target_data[j].is_separator()) {

                            }
                        }
                        break;
                    case stim::DEM_SHIFT_DETECTORS:
                        det_offset += instruction.target_data[0].data;
                        break;
                    case stim::DEM_DETECTOR:
                        break;
                    case stim::DEM_LOGICAL_OBSERVABLE:
                        break;
                    default:
                        break;
                }
            }
        }
    }
};

pm::Graph detector_error_model_to_matching_graph(stim::DetectorErrorModel& detector_error_model) {
    return pm::Graph(10);
}


