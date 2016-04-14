program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::m = 320, n = 200
  COOP_INT::i
  COOP_REAL::a(m, max(m, n)), w(max(m, n)), v(n, n), acopy(m, max(m, n)), wfill(m, n)
  call random_number(a)
!  call coop_print_matrix(a, m, n)
  acopy = a
  wfill = 0.
  call coop_svd_decompose(m, n, a, w, v)
!  call coop_svdcmp(m, n, a, m, max(m, n), w, max(m, n), v, n, n)
  do i=1, min(m, n)
     wfill(i,i) = w(i)
  enddo
  a(1:m, 1:n) =  matmul(matmul(a(1:m, 1:m), wfill), transpose(v)) -  acopy(1:m, 1:n)
  print*
!  call coop_print_matrix(a, m, n)
  print*
  print*, sum(a(1:m, 1:n)**2)


end program Test  
