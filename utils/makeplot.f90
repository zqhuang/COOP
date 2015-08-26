program TestNpeak
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT,parameter::n=1024
  COOP_REAL x(0:n-1), y(0:n-1)
  COOP_COMPLEX c(0:n/2)
  COOP_INT i
  c(0:1) = 0.d0
  do i=2, n/2
     c(i) = coop_random_complex_gaussian()/i
  enddo
  call coop_set_uniform(n, x, 0.d0, 1.d0)
  call fft_1d_backward(n, c, y)
  do i=0, n-1
     write(*,*) x(i), y(i)
  enddo
  
  
end program TestNpeak
