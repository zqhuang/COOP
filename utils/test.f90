program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 200, m = 8
  COOP_INT i
  COOP_REAL x(n), y(n), yy(n), c(m), a, b
  a = 0.d0
  b = 8.d0
  call coop_set_uniform(n, x, a, b)
  y=cos(x)+x**2/15.d0
  call coop_chebfit(n, x, y, m, a, b, c)
  do i=1, n
     call coop_chebeval(m, a, b, c, x(i), yy(i))
     print*, x(i), y(i), yy(i)
  enddo
  
end program Test
