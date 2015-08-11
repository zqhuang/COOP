program Test
  use coop_wrapper_typedef
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  integer,parameter::n = 21
  integer,parameter::fit_order = 5
  real*8 xa(n), ya(n), y, dy, c(2*fit_order+1), dev
  type(coop_smooth_fit)::sf
  xa = (/ -1.d0, -0.8d0, -0.5d0, -0.2d0, -0.1d0, 0.d0, 0.1d0, 0.2d0, 0.28d0,0.31d0, 0.32d0, 0.4d0, 0.41d0, 0.5d0, 0.55d0, 0.65d0, 0.8d0, 0.85d0, 0.9d0, 1.2d0, 1.5d0 /)

  !!example of extrapolation of noisy data
  call random_number(ya)
  ya = sin(xa) +(ya-0.5d0)/1000.d0

  call sf%fit(n, xa, ya, 6)
  print*, sf%eval(2.d0), sin(2.d0)
  
  
end program Test
