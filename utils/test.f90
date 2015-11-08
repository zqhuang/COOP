program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 1024
  COOP_REAL::x(n), y(n)
  COOP_INT::i
  call coop_set_uniform(n, x, -2.d0, 8.d0)
  y = erfc(x)
  do i=1, n
     write(*,*) x(i), coop_InverseErfc(y(i)) - x(i)
  enddo
end program Test
