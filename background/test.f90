program test
  use coop_wrapper_typedef
  use coop_background
  implicit none
#include "constants.h"
  type(coop_cosmology_background) bg
  bg = coop_cosmology_background(h=COOP_REAL_OF(0.69d0))
  call bg%add_species(coop_baryon(COOP_REAL_OF(0.045d0)))
  call bg%add_species(coop_cdm(COOP_REAL_OF(0.26d0)))
  call bg%add_species(coop_radiation(COOP_REAL_OF(0.000042d0)))
  call bg%add_species(coop_neutrinos_massive(COOP_REAL_OF(0.001d0), COOP_REAL_OF(0.00003d0)))
  call bg%add_species(coop_neutrinos_massless(COOP_REAL_OF(0.00001d0)))
  call bg%add_species(coop_de_quintessence(bg%Omega_k, COOP_REAL_OF(0.3d0), COOP_REAL_OF(0.d0), COOP_REAL_OF(0.5d0)))
  call bg%print


end program test
