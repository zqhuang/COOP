program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::n=44
  type(coop_file)::fp
  COOP_REAL::cov(n, n), L(n, n), e(n)
  COOP_INT i
  call coop_MPI_init()
  call fp%open("testmat.dat","r")
  call coop_read_matrix(fp%unit, cov, n, n)
  print*, sum(abs(cov-transpose(cov))) !!test symmetry
  L = cov
  call coop_matsym_diagonalize(L, e, .true.)
  do i=1, n
     print*, i, e(i)
  enddo
  call fp%close()
  call coop_MPI_finalize()
  
end program TestNpeak
