program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_INT, parameter::n=200
  COOP_REAL x(n), y(n)
  type(coop_function)::f
  call coop_set_uniform(n, x, -2.d0, 0.d0)
  x = exp(x)
  y = x**2
  call f%init(n = n, xmin = x(1), xmax=x(n), f = y, xlog = .true., check_boundary = .false.)
  print*, f%eval(0.3d0), 0.3d0**2
end program Test
