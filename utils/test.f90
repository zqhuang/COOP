program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_covmat)::covmat
  COOP_INT, parameter::n = 30
  COOP_INT i
  COOP_REAL::converge_R, lambda
  call coop_MPI_init()
  call covmat%alloc(n)
  covmat%mean = 0.d0 !coop_MPI_Rank()/10.d0
  covmat%mult = 1.d0
  call random_number(covmat%c)
  covmat%c = 2.d0*(covmat%c + transpose(covmat%c))-2.d0
  do i=1, n
     covmat%c(i, i) = covmat%c(i, i) + 30.d0
  enddo
  call covmat%normalize()  
  call covmat%MPI_Sync(converge_R = converge_R)
  if(coop_MPI_rank().eq.0)then
     write(*,*) converge_R
  endif
  call coop_MPI_finalize()
  
end program TestNpeak
