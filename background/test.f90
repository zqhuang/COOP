program test
  use coop_wrapper_utils
  use coop_background_mod
  implicit none
#include "constants.h"
  type(coop_cosmology_background) bg
  type(coop_species) de
  call bg%init(h=COOP_REAL_OF(0.682d0))
  call bg%add_species(coop_baryon(COOP_REAL_OF(0.047d0)))
  call bg%add_species(coop_cdm(COOP_REAL_OF(0.253d0)))
  call bg%add_species(coop_radiation(bg%Omega_radiation()))
  call bg%add_species(coop_neutrinos_massive( bg%Omega_nu_from_mnu_eV(0.06d0),bg%Omega_massless_neutrinos_per_species()))
  call bg%add_species(coop_neutrinos_massless(bg%Omega_massless_neutrinos_per_species()*2.))
  de = coop_de_quintessence(bg%Omega_k(), 0.d0, 0.d0, 0.d0)
  call bg%add_species(de)
  call bg%setup_background()
  print*, de%wofa(0.d0), de%wofa(1.d0)
end program test
