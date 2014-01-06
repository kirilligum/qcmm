#include "dfdv_num.hpp"
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


dfdv_num::dfdv_num(
  std::function<double(std::vector<double>)> f,
  double h) 
  : f(f), h(h), div_2h(0.5/h) {
  }

std::vector<double> dfdv_num::operator()(std::vector<double> &v) {
  std::vector<double> fv (v.size(),0.0);
  boost::for_each(v,fv,[&v,this](double &x, double &fi) {
      auto tmp = x;
      x+=h;
      auto fp =f(v);
      x = tmp -h;
      auto fm =f(v);
      x = tmp;
      fi =  (fp-fm)*div_2h;
      });
  return fv;
}

std::vector<double> dfdv_num::operator()(const std::vector<double> &cv) {
  std::vector<double> v(cv);
  return operator()(v);
}
