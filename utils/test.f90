program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL,dimension(:),allocatable::a
  call coop_MPI_init()
  allocate(a(2))
  a = 1.d0
  call coop_MPI_sum(a)
  write(*,*) coop_MPI_rank(), a
  
  call coop_MPI_finalize()
  
end program TestNpeak
