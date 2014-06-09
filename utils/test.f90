program plot
  use coop_wrapper_utils
  use coop_ode_mod
  use coop_random_mod
  implicit none
#include "constants.h"
  integer, parameter::n=5
  real x(n)
  real indices(n)
  integer i, j
  call coop_random_init()
  call random_number(x)
  indices = (/ (i, i=1,n) /)
  do i=1, n
     print*, i, indices(i), x(i)
  enddo
  print*
  call coop_quicksortacc(x, indices)
  do i=1, n
     print*, i, indices(i), x(i)
  enddo
end program plot
