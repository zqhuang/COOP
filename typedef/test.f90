program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_REAL::eta = -0.4d0, t0 = 0.9d0
  COOP_INT,parameter::n=128
  COOP_REAL::weff(n), z(n), a(n)
  COOP_INT::i
  call coop_set_uniform(n, z, 0.d0, 30.d0)
  a = 1.d0/(1.d0+z)
  do i=1, n
     weff(i) = page_effective_wdeofa(t0, eta, a(i))
     write(*,*)  z(i), weff(i)     
  enddo
end program Test  
