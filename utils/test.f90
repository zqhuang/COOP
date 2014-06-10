program plot
  use coop_wrapper_utils
  use coop_ode_mod
  use coop_random_mod
  use coop_sort_mod
  implicit none
#include "constants.h"
  integer, parameter::n= 31
  integer, parameter::nbins=50
  integer,allocatable:: ind(:)
  real*8,allocatable::x(:)
  real center(nbins), density(nbins)
  integer bins_used, i
  call coop_random_init()
  allocate(x(n), ind(n))
  call random_number(x)
  x(100) = -2.
  x(1001) = 2.

  call coop_quicksort(x)
  do i=1, n, max(n/100,1)
     print*, i, x(i)
  enddo
  
  call coop_quicksort(x)
  do i=1, n, max(n/100,1)
     print*, i, x(i)
  enddo

end program plot
