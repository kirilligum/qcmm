#include "eom.hpp"

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

void eom::operator () ( const std::vector<double> &v, std::vector<double> &dvdt, const double ) {
  auto dfdv_at_v = dfdv(v);
  dvdt[0] = 0.0; ///> the first element is time that will be calculated later
  for (size_t i = 1; i < dfdv_at_v.size(); ++(++i)) {
    dvdt[i]   =   dfdv_at_v[i+1];
    dvdt[i+1] = - dfdv_at_v[i]  ;
  }
}

