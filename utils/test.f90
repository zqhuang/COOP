program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_list_realarr) rl
  integer,parameter::m  = 3, n = 10
  integer i
  real  x(m)
  do i= 1, n
     call random_number(x)
     call rl%push(x)
  enddo
  do i = 1, rl%n
     print*, rl%element(i)
  enddo
  call rl%sort(3)
  print*
  do i = 1, rl%n
     print*, rl%element(i)
  enddo


end program Test
