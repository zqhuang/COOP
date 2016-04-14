program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_INT, parameter::m = 22, n = 15
  COOP_INT::i
  COOP_REAL::a(m, max(m, n)), w(max(m, n)), v(n, n), acopy(m, max(m, n)), wfill(m, n)
  call random_seed()
  call random_number(a)
  acopy = a
  wfill = 0.d0
  call coop_svd_decompose(m, n, a, w, v)
  print*, w
  do i=1, min(m, n)
     wfill(i,i) = w(i)
  enddo
  a(1:m, 1:n) =  matmul(matmul(a(1:m, 1:m), wfill), transpose(v)) -  acopy(1:m, 1:n)
  print*
  print*
  print*, sum(a(1:m, 1:n)**2)


end program Test  
