#include <random>
#include "rand_gen.h"

std::mt19937 & global_urng( ) {
static std::mt19937 u{}; return u;
}

void randomize( ) {
static std::random_device rd{};
global_urng().seed( rd() );
}

void set_seed( unsigned s) { 
    global_urng().seed(s); 
}

double rand_float(double from, double to) {
    static std::uniform_real_distribution<> d{};
    using parm_t = decltype(d)::param_type;
    return d( global_urng(), parm_t{from, to} ); 
}