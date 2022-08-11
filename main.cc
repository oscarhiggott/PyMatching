#include "namespaced_main.h"
#include <stdexcept>
#include <iostream>


int main(int argc, const char** argv) {
    try {
        return pm::main(argc, argv);
    } catch (std::invalid_argument &ex) {
        std::cerr << ex.what() << "\n";
        return EXIT_FAILURE;
    }
}
