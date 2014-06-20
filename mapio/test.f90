program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  character(LEN=*),parameter::mapdir = "act/"
  character(LEN=*),parameter::fitsfile = mapdir//"I60.fits"
  type(coop_fits_image_cea)::cf
  call cf%open(fitsfile)
  call cf%header%print
  call cf%write("copy.fits")
  call cf%free
  call cf%open("copy.fits")
  call cf%header%print
end program test
