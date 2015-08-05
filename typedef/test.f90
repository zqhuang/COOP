program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_INT,parameter::n = 8192
  COOP_REAL:: x(n)
  do i=1, n
     x(i) = coop_random_Gaussian()
  enddo
end program Test
