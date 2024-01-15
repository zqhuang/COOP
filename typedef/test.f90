program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"

  real*8::eta, t0, q0, omegam, omegal, h, omegar, astar
  omegam = 0.3158
  omegar = 0. !4.187e-5/h**2
  omegal = 1.-omegam - omegar
  t0 = lcdm_tofa(omegam = omegam, omegak = 0.d0, a = 1.d0)
  q0 = omegam/2. - omegal
  eta  = 1. - (1.+q0)/2.*3.*t0**2
  print*, t0, eta
end program Test  
