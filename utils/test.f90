program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL::a(3,3), R(3,3)
  a = 0.d0
  a(1,1) = 2.d0
  a(2,2) = 3.d0
  a(1,2) = -1.5d0
  a(2,1) = a(1,2)
  a(3,3) = 1.2d0
  call coop_matsym_diag(3,3,a,R)
  a = matmul(transpose(R), matmul(a, R))
  print*, a(:,1)
  print*, a(:,2)  
  print*, a(:,3)
end program Test
