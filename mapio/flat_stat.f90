program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_fits_image_cea)::map
  COOP_STRING::fname
  COOP_REAL::mean, upper, lower, tail
  call coop_MPI_Init()
  fname = trim(coop_InputArgs(1))
  if(trim(fname).eq."")then
     write(*,*) "Syntax: ./FStat FILENAME"
     stop
  endif
  call map%open(fname)
  mean = sum(map%image)/map%npix
  write(*,*) "size: "//COOP_STR_OF(map%nside(1))//" x "//COOP_STR_OF(map%nside(2))
  write(*,*) "mean = ", mean
  write(*,*) "rms = ", sqrt(sum((map%image-mean)**2/map%npix))
  tail = 0.1585
  call array_get_threshold_double(map%image, map%npix, 1.-tail, lower)
  call array_get_threshold_double(map%image, map%npix, tail, upper)
  write(*,*) "1sigma lower, upper = ", lower, upper
  tail = 0.023
  call array_get_threshold_double(map%image, map%npix, 1.-tail, lower)
  call array_get_threshold_double(map%image, map%npix, tail, upper)
  write(*,*) "2sigma lower, upper = ", lower, upper
  tail = 0.0015
  call array_get_threshold_double(map%image, map%npix, 1.-tail, lower)
  call array_get_threshold_double(map%image, map%npix, tail, upper)
  write(*,*) "3sigma lower, upper = ", lower, upper
  write(*,*) "min max = ",  minval(map%image), maxval(map%image)
  write(*,*) "zero-value pixels: "//trim(coop_num2str(100.*count(map%image .eq. 0.d0)/dble(map%npix),"(F10.3)"))//"%"
end program test
