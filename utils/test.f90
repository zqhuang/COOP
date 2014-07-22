program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  integer,parameter::n = 200
  COOP_REAL x(n), y(n), k, b, r
  integer i
  call coop_set_uniform(n, x, 10.d0, coop_pi)
  call random_number(y)
  y = (2.*y-1.)
  y = x*2.+5.+y
  call coop_linear_least_square_fit(n, x, y, k, b, r)
  print*, k, b, r
  
end program Test
