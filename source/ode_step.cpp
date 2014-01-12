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
  typedef runge_kutta_cash_karp54<std::vector<double>> error_stopper_type;
  typedef controlled_runge_kutta<error_stopper_type> controlled_stepper_type;
  //controlled_stepper_type controlled_stepper(default_error_checker<double>(1.0e-6,1.0e-6,1.0,1.0));
  controlled_stepper_type controlled_stepper(default_error_checker<double>(1.0e-9,1.0e-9,1.0,1.0));
  //integrate_n_steps(controlled_stepper,eom(d_ham),v,0.0,time_step,(int)(time_duration/time_step),streamer(result));
  integrate_adaptive(controlled_stepper, eom, current, a[0],a[0]+dt,dt*1e-3);
  //integrate(eom, current, a[0],a[0]+dt,dt*1e-3);
  current[0]+=dt;
  return current;
}
