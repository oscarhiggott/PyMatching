#ifndef PYMATCHING2_STIM_IO_H
#define PYMATCHING2_STIM_IO_H

#include "stim.h"
#include "graph.h"


//stim::DetectorErrorModel::from_file();
//auto reader = stim::MeasureRecordReader::make(nullptr, /*stim::format_name_to_enum_map*/ stim::SAMPLE_FORMAT_B8, 0 ,0, 0);
//stim::SparseShot sparse_shot;
//reader->start_and_read_entire_record(sparse_shot);
////        stim::find_argument()

pm::Graph detector_error_model_to_matching_graph(stim::DetectorErrorModel& detector_error_model);




#endif //PYMATCHING2_STIM_IO_H
