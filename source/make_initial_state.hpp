#pragma once

#include <vector>
#include <random>

#include "ohmic_bath_param.hpp"
#include "pcet_param.hpp"


std::vector<double> make_initial_state(pcet_param sp, std::default_random_engine &gen) ;
