program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  
  integer,parameter::n=1024
  type(coop_function) f
  COOP_REAL x(n), y(n), xm
  logical,parameter::xlog = .true., ylog = .true.
  call coop_set_uniform(n, x, 0.001d0, 5.d0, logscale = xlog)     
  call random_number(y)
  y = exp(-x**2/2.d0)
  call f%init(n=n, xmin = x(1), xmax=x(n), f = y, xlog = xlog, ylog = ylog, check_boundary = .false., method = COOP_INTERPOLATE_QUADRATIC)
  call f%monotonic_solution(exp(-4.d0/2.d0), xm)
  print*, xm
end program Test
