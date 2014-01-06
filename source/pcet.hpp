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

typedef std::vector<double> vd; typedef std::vector<std::vector<double>> vvd; typedef std::vector<std::vector<std::vector<double>>> vvvd;

struct pcet {
  double e1,e2,v12,w1,w2,mr;
  electronic_aa_cart e;
  double operator()(const vd& v) ;
  void operator () ( const std::vector<double> &v, std::vector<double> &dvdt, const double ) ;
};
