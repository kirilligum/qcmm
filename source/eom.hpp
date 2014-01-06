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
//#include <boost/range/algorithm.hpp>
//#include <boost/range/algorithm_ext.hpp>
//#include <boost/range/numeric.hpp>
//#include <boost/math/constants/constants.hpp>

struct eom {
  std::function<std::vector<double>(std::vector<double>)> dfdv;
  void operator () ( const std::vector<double> &v, std::vector<double> &dvdt, const double ) ;
};

