#include <iostream>
#include <stdexcept>

#include "pymatching/fill_match/driver/namespaced_main.h"

int main(int argc, const char** argv) {
    try {
        return pm::main(argc, argv);
    } catch (std::invalid_argument& ex) {
        std::cerr << ex.what() << "\n";
        return EXIT_FAILURE;
    }
}
