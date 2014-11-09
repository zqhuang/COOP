program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_INT, parameter::n=20000
  COOP_SINGLE x(0:n), y(0:n)
  COOP_INT i
  type(coop_function)::f
  call coop_set_uniform(n+1, x, 0.5, 3.)
  y = sin(1000.d0*x)**2/x
  call coop_smooth_data(n+1, y, 20)
  do i=0, n
     print*, x(i), y(i)*x(i)
  enddo
end program Test
