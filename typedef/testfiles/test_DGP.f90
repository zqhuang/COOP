program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  real*8:: omegam, omegak, w0, wa, h, omegal, rho, p
  integer,parameter::n=300
  real*8::eta, t0, q0
  integer i
  real*8::a(n), z(n), dl(n), dla(n), err(n)
  omegak = 0.d0
  omegam = 0.3d0
  h = 0.7d0
  !!---------------------
  t0 = DGP_tofa(omegam, 1.d0)
  q0 = -0.3076923
  eta = 1.d0 - (1.d0+q0)*(1.5d0*t0**2)
  print*, t0, eta
  call coop_set_uniform(n, z, 0.01d0, 2.5d0)
  a=1.d0/(1.d0+z)
  do i=1, n
     dl(i) = DGP_dlofa(omegam, omegak, a(i))
     dla(i) = page_dlofa(t0, eta, omegak, a(i))
     err(i) = abs((dla(i)-dl(i))/dl(i))
  enddo
  print*,maxval(err)
end program Test  
