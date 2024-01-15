program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  integer,parameter::dof=8
  integer,parameter::n = dof*20  
  real*8 xa(n), ya(n), x, r(n)
  type(coop_smooth_fit)::sf
  COOP_INT::i

  
  call coop_set_uniform(n, xa, 0.d0, 1.d0)

  print*, "extrapolating smooth data"  
  do i=1, n
     ya(i) = f(xa(i))
  enddo
  call sf%fit(n, xa, ya, dof = dof)
  do i = 1, 20
     x=1.d0+0.1d0*i
     print*, x, sf%eval(x), f(x)
  enddo

  print*, "extrapolating noisy data"
  
  call random_number(r)
  ya = ya + (r-0.5d0)*0.001d0  !!add some noise
  call sf%fit(n, xa, ya, dof = dof)
  do i = 1, 20
     x=1.d0+0.1d0*i
     print*, x, sf%eval(x), f(x)
  enddo

  
contains

  function f(x)
    COOP_REAL::f, x
    f = x*sqrt(x**2+1.d0)*cos(x)
  end function f
  
end program Test
