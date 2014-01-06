#include "ode_step.hpp"
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
#include <boost/range/numeric.hpp>
//#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>

ode_step::ode_step(double dt_,std::function<void(const std::vector<double>&,std::vector<double>&,const double)> eom_): dt(dt_), eom(eom_){}

std::vector<double> ode_step::operator()(std::vector<double> a, std::vector<double> b) {
  using namespace boost::numeric::odeint;
  std::vector<double> current(a);
  integrate(eom, current, a[0],a[0]+dt,dt*1e-3);
  current[0]+=dt;
  return current;
}
