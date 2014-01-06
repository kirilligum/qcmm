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

struct ode_step {
  const double dt;
  std::function<void(const std::vector<double>&,std::vector<double>&,const double)> eom;
  ode_step(double dt_,std::function<void(const std::vector<double>&,std::vector<double>&,const double)> eom_);
  std::vector<double> operator()(std::vector<double> a, std::vector<double> b);
};
