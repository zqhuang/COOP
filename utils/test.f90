program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_asy)::fig
  COOP_INT, parameter::n=100
  COOP_INT i
  COOP_REAL::x(n), y(n),yapp(n)
  call coop_set_uniform(n, x, 0.d0, coop_SI_degree*4.d0*50.d0)
  yapp = 1.d0-x**2/4.d0+x**4/64.d0
  do i=1, n
     y(i) = coop_BessJ0(x(i))
     write(*,*) x(i), y(i), yapp(i)
  enddo

end program Test
