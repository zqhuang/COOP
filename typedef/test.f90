program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  
  integer,parameter::n=8192
  type(coop_function) f
  COOP_REAL x(n), y(n), xm
  call coop_set_uniform(n, x, log(0.1d0), log(10.d0))
  x = exp(x)
  y = x**2*exp(-x)
  call f%init(n=n, xmin = x(1), xmax=x(n), f = y, xlog = .true., ylog = .true., check_boundary = .false.)
  xm = f%maxloc()
  print*, f%derivative(2.d0)
  print*, xm, f%derivative(xm), f%eval(xm - 1.d-5) - f%eval(xm), f%eval(xm+1.d-5)-f%eval(xm)
end program Test
