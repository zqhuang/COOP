program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  character(LEN=*),parameter::mapdir = "act/"
  type(coop_fits_image_cea)::cf
  COOP_REAL::smooth_pixsize = 5. * coop_SI_arcmin
  call coop_random_init()
  call cf%open(mapdir//"I6.fits")
  call cf%get_flatmap(smooth_pixsize)
  call cf%simulate_flat(lmin=200, lmax=2500, Cls_file = "cls.txt")
  call cf%plot("simu/act");
end program test
