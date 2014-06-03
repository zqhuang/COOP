program test
  use coop_wrapper_typedef
  use coop_background_mod
  implicit none
#include "constants.h"
  type(coop_cosmology_background) bg
  bg = coop_cosmology_background(h=COOP_REAL_OF(0.72d0))
  call bg%add_species(coop_baryon(COOP_REAL_OF(0.045d0)))
  call bg%add_species(coop_cdm(COOP_REAL_OF(0.26d0)))
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massive(COOP_REAL_OF(0.001d0),bg%Omega_massless_neutrinos_per_species()*3))
  call bg%add_species(coop_de_quintessence(bg%Omega_k(), COOP_REAL_OF(0.3d0), COOP_REAL_OF(0.d0), COOP_REAL_OF(0.5d0)))
  call bg%print()
  print*, bg%Tnu()
  print*, bg%Omega_nu_from_mnu_eV(93./3.*0.05d0)*bg%h()**2
  
end program test
