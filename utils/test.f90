program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 128
  COOP_REAL::x(n), y(n)
  COOP_INT::i
  print*, coop_Gaussian_nu_of_P(0.8485d0), coop_Gaussian_nu_of_P(0.977d0)
  
  call coop_set_uniform(n, x, -5.d0, 5.d0)
  y = erf(x)
  do i=1, n
     write(*,*) x(i), coop_InverseErf(y(i)) - x(i), y(i)
  enddo
end program Test
