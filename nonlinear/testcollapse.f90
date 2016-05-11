program test
  use coop_wrapper_firstorder
  use coop_ellipse_collapse_mod
  implicit none
#include "constants.h"
  COOP_REAL::a(3), bprime(3)
  read(*,*) a
  call coop_ellipse_collapse_compute_bprime(a, bprime)
  print*, bprime
end program test
