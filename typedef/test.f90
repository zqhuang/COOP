program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  
  integer,parameter::n=1024
  type(coop_function) f
  COOP_REAL x(n), y(n), xm
  logical,parameter::xlog = .false., ylog = .false.
  call coop_set_uniform(n, x, 0.d0, 20.d0, logscale = xlog)     
  call random_number(y)
  y = tanh(x-5.d0)+2.d0
  call f%init(n=n, xmin = x(1), xmax=x(n), f = y, xlog = xlog, ylog = ylog, check_boundary = .false., method = COOP_INTERPOLATE_LINEAR)
  xm = 1.d0
  print*, f%method
  do while(xm .gt. 0.d0)
     print*, f%derivative(xm),  1.d0/cosh(xm-5.d0)**2
     print*, f%derivative2(xm), -2.d0*sinh(xm-5.d0)/cosh(xm-5.d0)**3
     read(*,*) xm
  enddo
end program Test
