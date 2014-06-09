program plot
  use coop_wrapper_utils
  use coop_ode_mod
  implicit none
#include "constants.h"
  integer, parameter::n=5
  doubleprecision x(n, n)
  integer i, j
  call random_number(x)
  call quicksort_double(x, n**2)
  do i=1, n
     do j=1, n
        print*, i, j,  x(j, i)
     enddo
  enddo
end program plot
