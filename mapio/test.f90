program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  implicit none
#include "constants.h"

  character(LEN=*),parameter::mapdir = "../data/cmb/maps/act/"
  character(LEN=*),parameter::fitsfile = mapdir//"Q50.fits"
  type(coop_fits_image_cea)::cf
  COOP_REAL output(2)
  COOP_INT input(2)

  call cf%open(fitsfile)
  call cf%header%print()
  print*, cf%n
  call cf%get_data()
end program test
