program testcosmo
  use coop_wrapper
  implicit none

#include "constants.h"

  COOP_COSMO_PARAMS = coop_arguments ( r= (/ 0.022032d0, 0.12038d0, 0.01041190d0, 0.0925d0, 0.d0, 0.06d0, 0.d0, 0.d0, 0.0d0, 0.d0 /) , i=(/ COOP_DE_COSMOLOGICAL_CONSTANT, 7, 4 /) )
  call coop_setup_global_cosmology()
  print*, coop_global_cosmology_cosmomc_theta(), COOP_COSMO%h()
  
end program testcosmo
