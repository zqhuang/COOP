program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  integer i
  COOP_REAL::t0
  do i= -10, 10
     t0 = 1. + 0.02*i
     print*, 1.-1.5*(t0)**2, t0
  enddo
end program Test
