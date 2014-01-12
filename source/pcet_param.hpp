#pragma once

#include "units_const.hpp"

struct pcet_param {
  units_const u;
  double eps_r_ev = 1.0, w_r_wn = 3000, min_r_ang =  0.0,
         eps_k_ev = 0.0, w_k_wn = 3000, min_k_ang = -0.5,
                         w_0_wn = 3000, min_0_ang = -0.5;
  double delta=0.03, m= 1836.1,
         eps_r = u.electron_volts_to_AU*eps_r_ev, w_r = u.wavenumbers_to_AU*w_r_wn, min_r = u.angstroms_to_AU*min_r_ang, k_r = m*w_r*w_r, 
         eps_k = u.electron_volts_to_AU*eps_k_ev, w_k = u.wavenumbers_to_AU*w_k_wn, min_k = u.angstroms_to_AU*min_k_ang, k_k = m*w_k*w_k, 
                                                  w_0 = u.wavenumbers_to_AU*w_0_wn, min_0 = u.angstroms_to_AU*min_0_ang, k_0 = m*w_0*w_0,
         n_shift = 0.366, n1=1.0,n2=1.0-n1,qr_initial=-0.5,pr_initial=0.0;
};
