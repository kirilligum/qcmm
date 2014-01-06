#pragma once
#include <vector>
#include <functional>

struct dfdv_num {
  std::function<double(std::vector<double>)> f;
  double h;
  double div_2h;
  dfdv_num( std::function<double(std::vector<double>)> f, double h = 1e-1);
  std::vector<double> operator()(std::vector<double> &v) ;
  std::vector<double> operator()(const std::vector<double> &cv) ;
};
