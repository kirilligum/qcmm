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
    //std::cout <<"muahaha\n";
  }

std::vector<double> dfdv_num::operator()(std::vector<double> &v) {
  //std::cout << "lalala" << std::endl;
  std::vector<double> fv (v.size(),0.0);
  boost::for_each(v,fv,[&v,this](double &x, double &fi) {
      auto tmp = x;
      x+=h;
      //std::cout << "v = "; std::copy(std::begin(v), std::end(v),std::ostream_iterator<double>(std::cout," ")); std::cout << std::endl;
      auto fp =f(v);
      //std::cout << "fp = " << fp << std::endl;
      x = tmp -h;
      //std::cout << "v = "; std::copy(std::begin(v), std::end(v),std::ostream_iterator<double>(std::cout," ")); std::cout << std::endl;
      auto fm =f(v);
      //std::cout << "fm = " << fm << std::endl;
      x = tmp;
      //std::cout << "fp - fm = " << fp-fm << std::endl;
      //std::cout << "2h  = " << 2*h << "  " << div_2h << std::endl;
      fi =  (fp-fm)*div_2h;
      });
  //std::cout << "fv = "; std::copy(std::begin(fv), std::end(fv),std::ostream_iterator<double>(std::cout," ")); std::cout << std::endl;
  return fv;
}

std::vector<double> dfdv_num::operator()(const std::vector<double> &cv) {
  std::vector<double> v(cv);
  return operator()(v);
}
