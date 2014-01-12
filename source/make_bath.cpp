#include "make_bath.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <limits>
#include <tuple>
#include <iterator>
#include <random>
//#include <mpi.h>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>
//#include <boost/math/constants/constants.hpp>

#include "ohmic_bath_param.hpp"

double j_bath(double eta, double w_c, double w) {
  return eta*w*exp(-w/w_c);
}

std::tuple<
std::vector<double> ,
std::vector<double> ,
std::vector<double> 
> make_bath(const ohmic_bath_param& bp, std::default_random_engine &gen) {
  std::vector<double> bath;
  std::vector<double> w(bp.modes),c(bp.modes);
  double dw = bp.wmax/bp.modes;
  for (size_t i = 0; i < bp.modes; ++i) {
    w[i]=(i+0.5)*dw;
    c[i]=sqrt(2*dw*j_bath(bp.eta,bp.wc,w[i])*bp.m*w[i]/M_PI);
    double ai = 2*tanh(bp.beta*w[i]*0.5),
      mean_q = -c[i]/(w[i]*w[i])/bp.m,
      mean_p = 0.0,
      stddev_q = 1/sqrt( ai*w[i]*bp.m ),
      stddev_p = sqrt( w[i]*bp.m/ai );
    bath.emplace_back( std::normal_distribution<double>(mean_q, stddev_q )(gen));
    bath.emplace_back( std::normal_distribution<double>(mean_p, stddev_p )(gen));
  }
  return std::make_tuple(bath,w,c);
}
