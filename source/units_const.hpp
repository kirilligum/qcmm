#pragma once

struct units_const {
  const double	
      picoseconds_to_AU = 41356.49296939619520264,
      femptoseconds_to_AU = picoseconds_to_AU/1e3,
      wavenumbers_to_AU = 4.55633525e-06, // this is hc in Hartrees*cm
      angstroms_to_AU = 1/0.5291772192,
      electron_volts_to_AU = 1/27.211,
      boltzman_constant = 3.16681151E-06; // this is k in Hartrees/K
};
