program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  integer,parameter::n = 20
  COOP_REAL a(n, n), w(n), v(n, n), acopy(n, n)
  integer i
  call random_number(a)
  acopy = a
  call coop_matrix_sorted_svd(n, a, w, v)
  do i=1, n
     a(:,i)=a(:, i)*w(i)
  enddo

  a = acopy - matmul(a, transpose(v))
  print*, maxval(a), minval(a)

  print*, w
end program Test
