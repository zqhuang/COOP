program Test
  use coop_wrapper_typedef
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::n = 500000
  COOP_REAL,dimension(:),allocatable::x
  COOP_REAL::xbar, sigma, A
  COOP_INT::i
  call coop_MPI_init()
  call coop_random_init()
  allocate(x(n))
  do i=1, n
     x(i) = coop_random_Gaussian()*2.d0 - 1.5d0 + 0.02 * coop_random_Gaussian()**3
  enddo
  call coop_fit_Gaussian(x, n/500, xbar, sigma, A)
  print*, xbar, sigma, A
  call coop_asy_histogram(x, n/1000, xlabel = "$x$", ylabel = "$y$", filename = "x.txt", fit_gaussian = .true.)
  call coop_MPI_Finalize()

end program Test
