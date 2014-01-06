#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>

#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
using namespace std;

struct dfdv_num {
  function<double(vector<double>)> f;
  double h;
  double div_2h;
  dfdv_num(
    function<double(vector<double>)> f,
    double h = 1e-9) : f(f), h(h), div_2h(0.5/h) {}
  vector<double> operator()(vector<double> &v) {
    std::vector<double> fv (v.size(),0.0);
    boost::for_each(v,fv,[&v,this](double &x, double &fi) {
        auto tmp = x;
        x+=h;
        auto fp =f(v);
        x = tmp -h;
        auto fm =f(v);
        x = tmp;
        fi =  (fp-fm)/h*0.5;
        });
    return fv;
  }
  vector<double> operator()(const vector<double> &cv) {
    vector<double> v(cv);
    return (v);
  }
};

double x2(const vector<double>& cv) {
  vector<double> v(cv);
  for(auto &i : v) i*=i*i;
  return accumulate(begin(v)+1,end(v),v.front());
}

int main(int argc, char const *argv[]) {
  cout << " hi \n";
  vector<double> v {3.0};
  cout << v[0] << " " << v[1] << endl;
  cout << "dfdv  ( " << dfdv_num{x2}(v)[0] << " )"<< endl;
  for(auto i : dfdv_num{x2,0.1}(v)) cout << i << " . "; cout << endl;
  return 0;
}
