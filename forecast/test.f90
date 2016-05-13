program test
  use coop_wrapper_firstorder
  implicit none
#include "constants.h"
  COOP_REAL::r
  COOP_INT::ic
  do ic = 1, 3
     r = sqrt(-2.d0*log(coop_IncompleteGamma(0.5d0, dble(ic)**2/2.d0)/coop_sqrtpi))
     print*, ic, r**2
  enddo
end program test
