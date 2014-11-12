program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_INT, parameter::n=1000
  COOP_REAL x(n), y(n), xx(n*10), yy1(n*10), yy2(n*10)
  COOP_INT i
  type(coop_function)::f
  call coop_set_uniform(n, x, 0.5d0, 3.d0)
  call coop_set_uniform(n*10, xx, 0.5d0, 3.d0)
  y = sin(x)/x
  call f%init_NonUniform(x, y, smooth_boundary = .true.)
  yy1 = -sin(xx)/xx - 2.*cos(xx)/xx**2 + 2*sin(xx)/xx**3
  do i=1, n*10
     yy2(i) = f%derivative2(xx(i))
     print*, xx(i), yy1(i), yy2(i)
  enddo

end program Test
