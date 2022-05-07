program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  real*8:: omegam, omegak, w0, wa, h, omegal, rho, p
  integer,parameter::n=300
  real*8::eta, t0, q0
  integer i, j
  real*8::a(n), z(n), dl(n), dla(n), err(n)
  omegak = -0.2585d0
  omegam = 0.3d0
  w0 = -0.8d0
  wa = 0.d0
  h = 0.6931d0
  !!---------------------
  do i=-100, 100
     do j=-100, 100
        omegam = 0.328d0 + 0.0001d0*i
        w0 = -0.803d0 + 0.0001d0 * j
        omegal = 1.d0- omegam - omegak
        rho = (omegam + omegal)*3.d0
        p = omegal*w0*3.d0
        t0  = w0wa_tofa(omegam = omegam, omegak = omegak, w0 = w0, wa=wa, a=1.d0)
        q0 = (rho+3.d0*p)/6.d0
        eta = 1.d0 - (1.d0+q0)*(1.5d0*t0**2)
        if(abs(t0 - 0.9719) .lt. 0.0001 .and. abs(eta - 0.2814).lt. 0.0001)&
             print*, omegam, w0, t0, eta
     enddo
  enddo
end program Test  
