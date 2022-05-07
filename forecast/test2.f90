program test
  use coop_wrapper_firstorder
  use coop_page_mod
  implicit none
#include "constants.h"
  real*8,parameter::omegam = 0.3d0
  real*8,parameter::h = 0.7d0  
  real*8,parameter::omegak = 0.d0  
  real*8,parameter::t0 = (0.2629/omegam)**0.277
  real*8::eta =1.d0- omegam*t0**2*2.25 
  integer:: i, n
  real*8::t, a, z
  n  = 100
  do i=1, n
     z = dble(i)/n * 2.d0
     print*, z, page_distance_moduli(t0, eta, omegak, h, z, z ) - lcdm_distance_moduli(omegam, omegak, h, z, z)
  enddo

end program test
