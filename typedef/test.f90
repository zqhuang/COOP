program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  type(coop_function)::f
  COOP_INT,parameter::n = 2000
  logical::xlog, ylog
  COOP_REAL::x(n), y(n), k, r, b
  call coop_set_uniform(n, x, -2.d0, 3.d0)
  call random_number(y)
  y = 2.4d0*x - 1.3d0 + (y-0.5d0)/5.d0
  call coop_linear_regression(n, x, y, k, r, b)
  print*, k, r, b
end program Test  
