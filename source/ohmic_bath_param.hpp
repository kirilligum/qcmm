#pragma once
#include <cmath>

#include "units_const.hpp"

struct ohmic_bath_param {
  units_const u;
  std::size_t modes=10;
  double temperature_kelvin=300, wc_wn=600, bath_mass_au=1836.1, alpha=24, wmax_factor=7;
  double beta = 1/u.boltzman_constant/temperature_kelvin, eta=M_PI*0.5*alpha, wc=u.wavenumbers_to_AU*wc_wn, wmax=wc*wmax_factor, m = bath_mass_au;
};

