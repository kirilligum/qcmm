#pragma once

#include <vector>
#include <random>

#include "ohmic_bath_param.hpp"

double j_bath(double eta, double w_c, double w) ;

std::tuple<
std::vector<double> ,
std::vector<double> ,
std::vector<double> 
> make_bath(const ohmic_bath_param& bp, std::default_random_engine &gen);
