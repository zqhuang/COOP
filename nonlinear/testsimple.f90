program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::n = 250
  COOP_INT::i
  COOP_REAL::theta(n), a(n), r(n)
  call coop_set_uniform(n, theta, 0.02d0, coop_2pi-0.02)
  r = (1.d0-cos(theta))/(theta-sin(theta))**(2.d0/3.d0)*2.d0/6.d0**(2.d0/3.d0)
  a = (3.d0/4.d0)**(2.d0/3.d0)*(theta-sin(theta))**(2.d0/3.d0)/(5.d0/3.d0)/1.6865d0
  do i=1, n
     write(*,*) a(i), r(i)*a(i)
  enddo
end program test
