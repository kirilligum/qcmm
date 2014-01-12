#pragma once
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
//#include <mpi.h>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
//#include <boost/range/numeric.hpp>
//#include <boost/math/constants/constants.hpp>

#include "electronic_aa_cart.hpp"
#include "ohmic_bath_param.hpp"
#include "pcet_param.hpp"

typedef std::vector<double> vd; typedef std::vector<std::vector<double>> vvd; typedef std::vector<std::vector<std::vector<double>>> vvvd;

struct pcet {
  pcet_param sp;
  ohmic_bath_param bp;
  std::vector<double> bw;
  std::vector<double> bc;
  electronic_aa_cart e;
  double operator()(const vd& v) ;
  void operator () ( const std::vector<double> &v, std::vector<double> &dvdt, const double ) ;
};
