program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  integer, parameter::n = 31
  COOP_REAL r(n)
  integer i
  call coop_random_init()
  do i=1,n
     r(i) = coop_random_Gaussian()
  enddo
  print*, sum(r**2)

end program Test
