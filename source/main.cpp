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
//#include <boost/numeric/odeint.hpp>
typedef std::vector<double> vd; typedef std::vector<std::vector<double>> vvd; typedef std::vector<std::vector<std::vector<double>>> vvvd;

struct electronic_aa_cart {
  double n_shift;
  double x(double n, double q) { return sqrt(2*(n+n_shift))*cos(q); }
  double p(double n, double q) { return -sqrt(2*(n+n_shift))*sin(q); }
  double n(double x, double p) { return (x*x+p*p)*0.5 - n_shift;}
};

int main(int argc, char const *argv[]) {
  using namespace std;
  cout << "hi\n";
  double n_shift = 0.366, n1=1.0,n2=0.0,qr=-0.5,pr=0.0;
  std::random_device rd;
  std::uniform_real_distribution<double> ran_pi(0.0,2*M_PI);
  double q1 = ran_pi(rd);
  double q2 = ran_pi(rd);
  double xn1 = electronic_aa_cart{n_shift}.x(n1,q1);
  double pn1 = electronic_aa_cart{n_shift}.p(n1,q1);
  double xn2 = electronic_aa_cart{n_shift}.x(n2,q2);
  double pn2 = electronic_aa_cart{n_shift}.p(n2,q2);
  vd initial_state = {0.0,xn1,pn1,xn2,pn2,qr,pr}; ///> time + nn + r
  std::copy(std::begin(initial_state), std::end(initial_state),std::ostream_iterator<double>(cout," ")); cout << endl;
  cout << electronic_aa_cart{n_shift}.n(xn1,pn1) << endl;
  cout << electronic_aa_cart{n_shift}.n(xn2,pn2) << endl;
  return 0;
}
