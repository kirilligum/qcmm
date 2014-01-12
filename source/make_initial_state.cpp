#include "make_initial_state.hpp"

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

#include "electronic_aa_cart.hpp"
#include "ohmic_bath_param.hpp"
#include "pcet_param.hpp"
#include "make_bath.hpp"

//typedef std::vector<double> vd; typedef std::vector<std::vector<double>> vvd; typedef std::vector<std::vector<std::vector<double>>> vvvd;

std::vector<double> make_initial_state(pcet_param sp, std::default_random_engine &gen) {
  using namespace std;
  std::uniform_real_distribution<double> ran_pi(0.0,2*M_PI);
  double q1 = ran_pi(gen);
  double q2 = ran_pi(gen);
  electronic_aa_cart eac{sp.n_shift};
  double xn1 = eac.x(sp.n1,q1);
  double pn1 = eac.p(sp.n1,q1);
  double xn2 = eac.x(sp.n2,q2);
  double pn2 = eac.p(sp.n2,q2);
  vector<double> initial_state = {0.0,xn1,pn1,xn2,pn2,sp.qr_initial,sp.pr_initial}; ///> time + nn + r
  //std::copy(std::begin(initial_state), std::end(initial_state),std::ostream_iterator<double>(cout," ")); cout << endl;
  //std::copy(std::begin(initial_state), std::end(initial_state),std::ostream_iterator<double>(cout," ")); cout << endl;
  return initial_state;
}
