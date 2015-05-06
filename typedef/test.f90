program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  type(coop_cosmology_background)::bg
  call bg%init(h = 0.68d0)
  print*, bg%omega_radiation() + bg%omega_massless_neutrinos()
  
end program Test
