program test
  use coop_wrapper_firstorder
  use coop_page_mod
  implicit none
#include "constants.h"
  real*8::w = -0.5d0
  real*8,parameter::omegam = 0.3d0
  real*8,parameter::h = 0.7d0  
  real*8,parameter::omegak = 0.d0  
  real*8::t0 
  real*8::eta, q0, dt
  integer:: i, n
  real*8::t, a, z

  do i=-3, 9
     w = -1.+0.1*i
     t0 = wcdm_tofa(omegam, omegak, w, 1.d0)
     eta =1.d0 -  1.5 * (omegam/2 + (1.+3.*w)/2. * (1.-omegam-omegak) + 1.) * t0**2
     print*, omegam , w, t0, eta

  enddo
  stop
  n  = 20
  do i=1, n
     z = dble(i)/n * 1.5d0
     print*, z, page_distance_moduli(t0, eta, omegak, h, z, z ) - wcdm_distance_moduli(omegam, omegak, w, h, z, z)
  enddo

end program test
