program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_covmat)::covmat
  call coop_MPI_init()
  call covmat%alloc(2)
  covmat%mean = coop_MPI_Rank()
  covmat%c(1, 1) = 1.d0+coop_MPI_Rank()
  covmat%c(2, 2) = 1.d0+coop_MPI_Rank()
  covmat%c(2, 1) = 0.4d0
  call covmat%normalize()
  covmat%mult = 1.d0
  call covmat%MPI_Sync()
  write(*,*) coop_MPI_Rank(), covmat%sigma, covmat%c
  call coop_MPI_finalize()
  
end program TestNpeak
