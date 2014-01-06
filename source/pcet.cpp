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
  return pr*pr*0.5/mr+
    (e1+w1*w1*mr*qr*qr*0.5)*n1+
    (e2+w2*w2*mr*qr*qr*0.5)*n2+
    v12*(xn1*xn2+pn1*pn2);
    //2*v12*sqrt((n1+e.n_shift)*(n2+e.n_shift))*(xn1*xn2+pn1*pn2); ///> TODO
}

void pcet::operator () ( const std::vector<double> &v, std::vector<double> &dvdt, const double ) {
  double xn1 = v[1],pn1 = v[2],xn2 = v[3],pn2 = v[4],qr = v[5],pr = v[6];
  double n1=e.n(xn1,pn1),n2=e.n(xn2,pn2);
  double mr_qr_qr_05 = mr*qr*qr*0.5;
  double w1_w1_mr_qr_qr_05 = w1*w1*mr_qr_qr_05;
  double w2_w2_mr_qr_qr_05 = w2*w2*mr_qr_qr_05;
  dvdt[0]=0.0;
  dvdt[1]=(e1+w1_w1_mr_qr_qr_05)*v[2]+v12*v[4];
  dvdt[2]=(e1+w1_w1_mr_qr_qr_05)*v[1]+v12*v[3];
  dvdt[2]*=-1;
  dvdt[3]=(e2+w2_w2_mr_qr_qr_05)*v[4]+v12*v[2];
  dvdt[4]=(e2+w2_w2_mr_qr_qr_05)*v[3]+v12*v[1];
  dvdt[4]*=-1;
  dvdt[5]=v[6]/mr;
  dvdt[6]=w1*w1*mr*qr*n1 + w2*w2*mr*qr*n2;
  dvdt[6]*=-1;
}
