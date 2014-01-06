#include "electronic_aa_cart.hpp"

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

double electronic_aa_cart::x(double n, double q) { return sqrt(2*(n+n_shift))*cos(q); }
double electronic_aa_cart::p(double n, double q) { return -sqrt(2*(n+n_shift))*sin(q); }
double electronic_aa_cart::n(double x, double p) { return (x*x+p*p)*0.5 - n_shift;}

