program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  integer, parameter::n = 31, m = 1000
  COOP_REAL r(n), chisq(m)
  integer i, j
  call coop_random_init()
  do j=1, m
     do i=1,n
        r(i) = coop_random_Gaussian()
     enddo
     chisq(j) = sum(r**2)
  enddo
  print*, maxval(chisq)
end program Test
