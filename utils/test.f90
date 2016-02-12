program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT::i
  do i= 1, 5
     write(*,*) sqrt(-2.d0*log(coop_IncompleteGamma(0.5d0, dble(i)**2/2.d0)/sqrt(coop_pi)))
  enddo


contains

  function gaussian(x)
    COOP_REAL::x, gaussian
    gaussian = exp(-x**2/2.d0)
  end function gaussian
end program Test
