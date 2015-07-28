program TestNpeak
  use,intrinsic::iso_c_binding
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n=31, m = 20  
  COOP_REAL::image(-n:n, -n:n, m)
  COOP_COMPLEX::fimage(0:n, 0:2*n, m)
  COOP_INT::i 
  call random_number(image)
  !$omp parallel do
  do i=1, m
     call coop_fft_forward(2*n+1, 2*n+1, image(:,:,i), fimage(:, :, i))
  enddo
  !$omp end parallel do
  print*, fimage(0, 0,1)
  print*, fimage(0, 1,1), fimage(0, 2*n,1)
end program TestNpeak
