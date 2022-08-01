#include "graph.h"
#include "graph_flooder.h"
#include "graph_fill_region.h"


pm::GraphFlooder::GraphFlooder(pm::Graph &graph) : graph(std::move(graph)), time(0) {};

void pm::GraphFlooder::create_region(DetectorNode *node) {
    auto region = new GraphFillRegion();
    auto alt_tree_node =  new AltTreeNode(region);
    region->alt_tree_node = alt_tree_node;
    do_region_created_at_empty_detector_node(*region, *node);
}

void pm::GraphFlooder::do_region_created_at_empty_detector_node(GraphFillRegion &region,
                                                                DetectorNode &detector_node) {
    detector_node.reached_from_source = &detector_node;
    detector_node.distance_from_source = 0;
    detector_node.region_that_arrived = &region;
    region.shell_area.push_back(&detector_node);
    reschedule_events_at_detector_node(detector_node);
}

void pm::GraphFlooder::reschedule_events_at_detector_node(DetectorNode &detector_node) {
    detector_node.invalidate_involved_schedule_items();

}