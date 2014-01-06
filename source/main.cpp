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

#include "dfdv_num.hpp"
#include "eom.hpp"
#include "electronic_aa_cart.hpp"
#include "pcet.hpp"
#include "ode_step.hpp"

//typedef std::vector<double> vd; typedef std::vector<std::vector<double>> vvd; typedef std::vector<std::vector<std::vector<double>>> vvvd;

int main(int argc, char const *argv[]) {
  using namespace std;
  cout << "hi\n";
  double n_shift = 0.366, n1=1.0,n2=1.0-n1,qr=-0.5,pr=0.0;
  std::random_device rd;
  std::uniform_real_distribution<double> ran_pi(0.0,2*M_PI);
  double q1 = ran_pi(rd);
  double q2 = ran_pi(rd);
  electronic_aa_cart eac{n_shift};
  double xn1 = eac.x(n1,q1);
  double pn1 = eac.p(n1,q1);
  double xn2 = eac.x(n2,q2);
  double pn2 = eac.p(n2,q2);
  vector<double> initial_state = {0.0,xn1,pn1,xn2,pn2,qr,pr}; ///> time + nn + r
  std::copy(std::begin(initial_state), std::end(initial_state),std::ostream_iterator<double>(cout," ")); cout << endl;
  double e1=1.0,e2=0.0,v12=1.0e1,w1=1.0,w2=1.0,mr=1.0;
  size_t time_steps = 4e2;
  double end_time =2e0;
  double dt = end_time/time_steps;
  vector<vector<double>> traj(time_steps);
  traj[0]=initial_state;
  cout << "anal energy = " << pcet{e1,e2,v12,w1,w2,mr,eac}(initial_state) << endl;
  cout << "num energy = " << pcet{e1,e2,v12,w1,w2,mr,eac}(initial_state) << endl;
  vector<double> ader(initial_state);
  pcet{e1,e2,v12,w1,w2,mr,eac}(initial_state,ader,1e-6);
  cout << "anal deriv = ";
  std::copy(std::begin(ader), std::end(ader),std::ostream_iterator<double>(std::cout," ")); std::cout << std::endl;
  vector<double> nder(initial_state);
  eom{dfdv_num{pcet{e1,e2,v12,w1,w2,mr,eac}}}(initial_state,nder,1e-6);
  cout << "num deriv = ";
  std::copy(std::begin(nder), std::end(nder),std::ostream_iterator<double>(std::cout," ")); std::cout << std::endl;
  boost::partial_sum(traj, traj.begin(),ode_step(dt,pcet{e1,e2,v12,w1,w2,mr,eac}));
  ofstream ot("traj.txt"); for(auto i:traj) { std::copy(std::begin(i), std::end(i),std::ostream_iterator<double>(ot," ")); ot << std::endl;}
  ofstream on("trajnne.txt"); for(auto i:traj) { 
    on << i[0] << "  ";
    on << eac.n(i[1],i[2]) << "  ";
    on << eac.n(i[3],i[4]) << "  ";
    on << pcet{e1,e2,v12,w1,w2,mr,eac}(i) << "  ";
    on << std::endl;
  }
  vector<vector<double>> trajn(time_steps);
  trajn[0]=initial_state;
  boost::partial_sum(trajn, trajn.begin(),ode_step(dt,eom{dfdv_num{pcet{e1,e2,v12,w1,w2,mr,eac}}}));
  ofstream otn("trajn.txt"); for(auto i:trajn) { std::copy(std::begin(i), std::end(i),std::ostream_iterator<double>(otn," ")); otn << std::endl;}
  ofstream onn("trajnnen.txt"); for(auto i:trajn) { 
    onn << i[0] << "  ";
    onn << eac.n(i[1],i[2]) << "  ";
    onn << eac.n(i[3],i[4]) << "  ";
    onn << pcet{e1,e2,v12,w1,w2,mr,eac}(i) << "  ";
    onn << std::endl;
  }
  return 0;
}
