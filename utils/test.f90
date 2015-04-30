program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_random_cycl)::cycl
  COOP_INT :: i, next
  call cycl%init(10)
  do i = 1, 20
     print*, cycl%next()
  enddo
  
end program TestNpeak
