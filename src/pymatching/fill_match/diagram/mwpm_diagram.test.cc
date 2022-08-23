#include "pymatching/fill_match/diagram/mwpm_diagram.h"

#include <fstream>

#include "gtest/gtest.h"

#include "stim.h"
#include "pymatching/fill_match/driver/mwpm_decoding.h"

TEST(mwpm_diagram, example) {
    stim::CircuitGenParameters params(45, 53, "memory");
    params.before_measure_flip_probability = 0.1;
    params.after_reset_flip_probability = 0.1;
    params.after_clifford_depolarization = 0.1;
    stim::Circuit circuit = stim::generate_rep_code_circuit(params).circuit;
    stim::DetectorErrorModel dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, true, true, false, 0, false, false);
    auto coords = pm::dem_detector_coords(dem);
    std::mt19937_64 rng(0);
    auto shot = stim::detector_samples(circuit, 1, false, false, rng);

    auto mwpm = pm::detector_error_model_to_mwpm(dem, 1000);
    std::vector<size_t> detection_events;
    for (size_t k = 0; k < coords.size(); k++) {
        if (shot[k][0]) {
            detection_events.push_back(k);
        }
    }
    detection_events.clear();
    detection_events.push_back(600);
    detection_events.push_back(610);
    detection_events.push_back(660);
    detection_events.push_back(700);
    detection_events.push_back(820);

    mwpm.flooder.queue.cur_time = 0;
    for (auto& detection : detection_events) {
        mwpm.create_detection_event(&mwpm.flooder.graph.nodes[detection]);
    }

    size_t id = 0;
    auto write_state = [&](pm::MwpmEvent ev){
        std::ofstream f;
        std::string s = std::to_string(id);
        while (s.size() < 5) {
            s.insert(s.begin(), '0');
        }
        f.open("/home/craiggidney/tmp/tmp" + s + ".svg");
        pm::write_mwpm_svg_diagram(coords, mwpm, ev, f);
        f.close();
        id++;
    };
    write_state(pm::MwpmEvent::no_event());
    size_t draw_step = 1000;
    auto next_draw_time = draw_step;

    while (true) {
        auto flood_event = mwpm.flooder.queue.dequeue();
        while (mwpm.flooder.queue.cur_time > next_draw_time) {
            auto t = mwpm.flooder.queue.cur_time;
            mwpm.flooder.queue.cur_time = next_draw_time;
            write_state(pm::MwpmEvent::no_event());
            mwpm.flooder.queue.cur_time = t;
            next_draw_time += draw_step;
        }
        if (!mwpm.flooder.dequeue_decision(flood_event)) {
            continue;
        }
        if (flood_event.tentative_event_type == pm::NO_FLOOD_CHECK_EVENT) {
            break;
        }

        auto mwpm_event = mwpm.flooder.process_tentative_event_returning_mwpm_event(flood_event);
        if (mwpm_event.event_type != pm::NO_EVENT) {
            write_state(mwpm_event);
            mwpm.process_event(mwpm_event);
            next_draw_time = mwpm.flooder.queue.cur_time;
        }
    }
    write_state(pm::MwpmEvent::no_event());

    ASSERT_TRUE(true);
}
