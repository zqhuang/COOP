program plot
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  integer, parameter::ny = 20, nx=30
  COOP_REAL fx(ny* nx)
  integer i, j
  COOP_COMPLEX fk(ny/2+1, nx)
  do i=1, nx*ny
     fx(i) = sin(i*0.1+0.2/i)+i/100.
  enddo
  print*, fx(1:10:2)
  call coop_fft_forward(nx, ny, fx, fk)
  call coop_fft_backward(nx, ny, fk, fx)
  print*, fx(1:10:2)
end program plot
