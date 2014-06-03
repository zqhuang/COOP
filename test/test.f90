program test
  use coop_wrapper_background
  implicit none
#include "constants.h"
  type(coop_cosmology_background) bg
  COOP_REAL, parameter::Omega_b = 0.045
  COOP_REAL, parameter::Omega_c = 0.26
  COOP_REAL, parameter::Omega_r = 4.e-5
  COOP_REAL, parameter::Omega_nu_massless = 3.e-5
  COOP_REAL, parameter::Omega_nu_massiv = 0.01
  call bg%init()
  call bg%add_species(coop_baryon(Omega_b))
  call bg%add_species(coop_cdm(Omega_c))
  call bg%print()
end program test
