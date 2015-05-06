program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_covmat)::covmat
  call coop_MPI_init()
  call covmat%alloc(2)
  covmat%mean = coop_MPI_Rank()+1.d0
  covmat%c(1, 1) = 1.d0+coop_MPI_Rank()
  covmat%c(2, 2) = 2.d0+coop_MPI_Rank()
  covmat%c(2, 1) = 0.4d0
  call covmat%normalize()
  covmat%mult = 1.d0
  if(coop_MPI_rank().eq.0)write(*,*) coop_MPI_Rank(), covmat%mean, covmat%sigma  
  call covmat%MPI_Sync()
  if(coop_MPI_rank().eq.0)write(*,*) coop_MPI_Rank(), covmat%mean, covmat%sigma
  call coop_MPI_finalize()
  
end program TestNpeak
