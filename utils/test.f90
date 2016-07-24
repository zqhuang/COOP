program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 4
  COOP_REAL::a(n,n), b(n, n)
  COOP_INT::i
  call coop_random_init()
  call random_number(b)
  a = b + transpose(b)
  do i=1, n
     a(i,i) = a(i,i) + 0.5d0+1.d0/i
  enddo
  b = a
  do i=1, n
     a(i, i+1:) = 0.d0
  enddo
  call coop_sympos_inverse(n, n, a, i)
  if(i.eq. 0)then
     a = matmul(a, b)
     call coop_print_matrix(a, n, n)
  else
     print*, i
  endif
end program Test  
