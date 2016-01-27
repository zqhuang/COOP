program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_fits_image_cea)::map, map2
  COOP_STRING::fname, fname2
  COOP_REAL::mean, rms, mean2, rms2
  fname = trim(coop_InputArgs(1))
  if(trim(fname).eq."")then
     write(*,*) "Syntax: ./FStat FILENAME"
     stop
  endif
  call map%open(fname)
  call map%simple_stat(mean=mean, rms=rms)
  fname2 = trim(coop_InputArgs(2))
  if(trim(fname2).ne."")then
     call map2%open(fname2)
     call map2%simple_stat(mean=mean2, rms=rms2)
     write(*,*) "correlation between two maps:", sum((map%image-mean)*(map2%image-mean2))/sqrt(sum((map%image-mean)**2)*sum((map2%image-mean2)**2))
     call map2%free()
  endif
  call map%free()
end program test
