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
#include "ohmic_bath_param.hpp"
#include "pcet_param.hpp"
#include "make_bath.hpp"
#include "make_initial_state.hpp"

//typedef std::vector<double> vd; typedef std::vector<std::vector<double>> vvd; typedef std::vector<std::vector<std::vector<double>>> vvvd;

int main(int argc, char const *argv[]) {
  using namespace std;
  cout << "hi\n";
  //std::random_device rd;
  size_t world_rank=1, seed = 1;
  std::default_random_engine gen(world_rank+seed);
  ohmic_bath_param  bp; std::vector<double> initial_bath,bw,bc;
  tie(initial_bath,bw,bc) = make_bath(bp,gen);
  pcet_param sp;
  auto initial_state = make_initial_state(sp,gen);
  initial_state.insert(end(initial_state),begin(initial_bath),end(initial_bath));
  std::copy(std::begin(initial_state), std::end(initial_state),std::ostream_iterator<double>(cout," ")); cout << endl;
  electronic_aa_cart eac{sp.n_shift};
  size_t time_steps = 4e2;
  double end_time =2e3;
  double dt = end_time/time_steps;
  vector<vector<double>> traj(time_steps);
  traj[0]=initial_state;
  boost::partial_sum(traj, traj.begin(),ode_step(dt,pcet{sp,bp,bw,bc,eac}));
  ofstream ot("traj.txt"); for(auto i:traj) { std::copy(std::begin(i), std::end(i),std::ostream_iterator<double>(ot," ")); ot << std::endl;}
  ofstream on("trajnne.txt"); for(auto i:traj) { 
    on << setw(20) << fixed << i[0] << "  ";
    on << setw(20) << fixed << eac.n(i[1],i[2]) << "  ";
    on << setw(20) << fixed << eac.n(i[3],i[4]) << "  ";
    on << setw(20) << fixed << pcet{sp,bp,bw,bc,eac}(i) << "  ";
    on << std::endl;
  }
  return 0;
}
