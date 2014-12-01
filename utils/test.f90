program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 2000
  COOP_INT i
  COOP_REAL x(n), y(n), a, b, f(n)
  COOP_COMPLEX z(n), ez(n)
  type(coop_asy)::fig
  a = 0.d0
  b = coop_2pi
  call coop_set_uniform(n, x, a, b)
  z= 2.d0*cmplx(cos(x), sin(x))
  ez = exp(2.d0*cos(x))*cmplx( cos(2.d0*sin(x)), sin(2.d0*sin(x)))
  f = abs(z*ez)-abs(0.25d0)
  do i=1, n
     print*, z(i), f(i)
  enddo
  
end program Test
