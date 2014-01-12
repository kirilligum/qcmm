#include "pcet.hpp"
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

#include "electronic_aa_cart.hpp"

double pcet::operator()(const vd& v) {
  double xn1 = v[1],pn1 = v[2],xn2 = v[3],pn2 = v[4],qr = v[5],pr = v[6];
  double n1=e.n(xn1,pn1),n2=e.n(xn2,pn2);
  double h_bath=0.0;
  for (size_t i = 7; i < v.size(); ++(++i)) {
    size_t ib=(i-7)/2;
    double mw2=bp.m*pow(bw[ib],2);
    h_bath+= v[i+1]*v[i+1]*0.5/bp.m
      + mw2*0.5*pow((v[i]+n1*bc[ib]/mw2),2);
  }
  double h11 = sp.eps_r + sp.w_r*sp.w_r*sp.m*pow((qr-sp.qr_initial),2)*0.5;
  double h22 = sp.eps_k + sp.w_k*sp.w_k*sp.m*pow((qr-sp.qr_initial),2)*0.5;
  double h12 = sp.delta;
  double h_avg= (h11+h22)*0.5;
  double h_delta= h11-h22;
  return pr*pr*0.5/sp.m+
    h_avg+
    h_delta*(n1-n2)*0.5+
    sp.delta*(xn1*xn2+pn1*pn2) + h_bath;
    //2*delta*sqrt((n1+e.n_shift)*(n2+e.n_shift))*(xn1*xn2+pn1*pn2); ///> TODO
}

void pcet::operator () ( const std::vector<double> &v, std::vector<double> &dvdt, const double ) {
  double xn1 = v[1],pn1 = v[2],xn2 = v[3],pn2 = v[4],qr = v[5],pr = v[6];
  double n1=e.n(xn1,pn1),n2=e.n(xn2,pn2);
  double mr_qr_qr_05 = sp.m*qr*qr*0.5;
  double w1_w1_mr_qr_qr_05 = sp.w_r*sp.w_r*mr_qr_qr_05;
  double w2_w2_mr_qr_qr_05 = sp.w_k*sp.w_k*mr_qr_qr_05;

  double h11 = sp.eps_r + sp.w_r*sp.w_r*sp.m*pow((qr-sp.qr_initial),2)*0.5;
  double h22 = sp.eps_k + sp.w_k*sp.w_k*sp.m*pow((qr-sp.qr_initial),2)*0.5;
  double h12 = sp.delta;
  double h_avg= (h11+h22)*0.5;
  double h_delta= h11-h22;
  double dh11 = sp.w_r*sp.w_r*sp.m*(qr-sp.qr_initial);
  double dh22 = sp.w_k*sp.w_k*sp.m*(qr-sp.qr_initial);
  double dh_avg= (dh11+dh22)*0.5;
  double dh_delta= dh11-dh22;

  double dbc=0;
  for (size_t i = 7; i < v.size(); ++(++i)) {
    size_t ib=(i-7)/2;
    double mw2=bp.m*pow(bw[ib],2);
    dbc+=bc[ib]*(v[i]+n1*bc[ib]/mw2);
    dvdt[i]=v[i+1]/bp.m;
    dvdt[i+1]=-mw2*v[i]-n1*bc[ib];
  }

  dvdt[0]=0.0;
  dvdt[1]= h_delta/2 * pn1 + h12 * pn2 + dbc * pn1;
  dvdt[2]= h_delta/2 * xn1 + h12 * xn2 + dbc * xn1;
  dvdt[2]*=-1;
  dvdt[3]= -h_delta/2 * pn2 + h12 * pn1;
  dvdt[4]= -h_delta/2 * xn2 + h12 * xn1;
  dvdt[4]*=-1;
  dvdt[5]=v[6]/sp.m;
  dvdt[6]= dh_avg + dh_delta*0.5*(n1-n2);
  dvdt[6]*=-1;
}
