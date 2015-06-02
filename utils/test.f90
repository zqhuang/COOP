program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::n=300000
  COOP_REAL::x(n)
  call coop_MPI_init()
  call coop_random_init()
  call coop_random_get_gaussian_vector(n, x)
  call coop_asy_histogram(x = x, nbins = 20, filename="test.txt")
  call coop_MPI_finalize()
  
end program TestNpeak
