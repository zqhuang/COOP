program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"

  real*8::eta, t0, q0, omegam, omegal, h, omegar, astar
  astar = 1.d0/1090.d0
  h = 0.68019908
  t0 = 0.9543577255135702
  omegam = 0.31062024
  omegar = 4.187e-5/h**2
  omegal = 1.-omegam - omegar
  q0 = omegam/2. - omegal
  eta  = 1. - (1.+q0)/2.*3.*t0**2
  print*, t0, eta
  print*, page_chiofa(t0, eta, astar), lcdm_chiofa(omegam, 0.d0, astar)
end program Test  
