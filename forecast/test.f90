program test
  use coop_wrapper_firstorder
  use coop_page_mod
  implicit none
#include "constants.h"
  real*8,parameter::omegam = 0.3d0
  real*8,parameter::h = 0.7d0  
  real*8,parameter::omegak = 0.7d0  
  real*8::t0 
  real*8::eta, omegal, q0, dt
  integer:: i, n
  real*8::t, a, z
  t0 = lcdm_tofa(omegam, omegak, 1.d0)
  omegal = 1.d0-omegam - omegak
  q0 = omegam/2.d0 - omegal
  eta =1.d0 -  1.5 * (1.d0 + q0) * t0**2
  dt = 0.01d0/lcdm_Hofa(omegam, omegak, 0.995d0)/0.995d0
  t = page_tofa(t0, eta, 0.995d0)
  n  = 100
  do i=1, n
     z = dble(i)/n * 2.d0
     print*, z, page_distance_moduli(t0, eta, omegak, h, z, z ) - lcdm_distance_moduli(omegam, omegak, h, z, z)
  enddo

end program test
