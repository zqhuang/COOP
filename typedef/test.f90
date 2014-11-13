program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_INT, parameter::n=1000
  COOP_REAL x(n), y(n), xx(n*10), yy1(n*10), yy2(n*10)
  COOP_INT i
  type(coop_function)::f
  COOP_REAL, parameter::a=0.1d0, b = 2.d0
  call coop_set_uniform(n, x, a, b)
  call coop_set_uniform(n*10, xx, a, b)
  y = exp(1/x)
  call f%init_NonUniform(x, y, ylog=.true.)
  yy1 = (1.d0/xx**4 + 2.d0/xx**3)*exp(1.d0/xx)
  do i=1, n*10
     yy2(i) = f%derivative2(xx(i))
     print*, xx(i), yy1(i), yy2(i)
  enddo

end program Test
