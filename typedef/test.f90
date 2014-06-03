program test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_INT,parameter::n=301
  COOP_INT i
  COOP_REAL, parameter::xmax = 6.d0
  COOP_REAL,dimension(n)::x, nx, y, y2, yp, fp
  COOP_REAL dx
  call coop_set_uniform(n, x, 0.d0, xmax)
  y = cos(x)
  yp = -sin(x)
  dx = x(2)-x(1)
  nx = x/dx
  call coop_spline_uniform(n, y, y2)
  do i=1, n
     call coop_splint_derv_uniform(n, x(1), dx, y, y2, x(i), fp(i))
  enddo
  do i=1, n
     write(*,*) x(i), fp(i), yp(i)
  enddo
end program test
