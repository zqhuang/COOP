program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::m = 100
  COOP_INT,parameter::n=40 !45+25
  type(coop_file)::fp
  COOP_REAL::cov(n, n), L(m, n), e(n)
  COOP_INT i, info
  call coop_MPI_init()
  call coop_import_matrix("testmat.dat", cov, n, n)
  L = 1.d0
  L(1:n, 1:n) = cov
  call coop_matsym_diagonalize(m, n, L, e)
  do i=1, n
     print*, i, e(i)
  enddo
  call coop_MPI_finalize()
  
end program TestNpeak
