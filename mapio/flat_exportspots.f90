program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  integer,parameter::lmin = 200
  integer,parameter::lmax = 2500
  character(LEN=*),parameter::mapdir = "act/"
  character(LEN=*),parameter::fitsfile = mapdir//"I6.fits"
  character(LEN=*),parameter::maskfile = mapdir//"WI6.fits"
  type(coop_fits_image_cea)::cf, mask, src
  type(coop_asy)::asy
  integer, parameter::n=300
  integer ix, iy, i, l
  type(coop_file) fp
  COOP_REAL map(n, n), Cls(lmin:lmax)
  COOP_REAL, parameter::smooth_scale = coop_SI_arcmin * 1.5

  call cf%open(fitsfile)
  call mask%open(maskfile)
  cf%image = (cf%image)*mask%image/maxval(mask%image)
  call cf%regularize(0.003d0)
  call cf%get_flatmap(smooth_scale)
  call cf%smooth_flat(lmin = lmin, lmax = lmax)
  call cf%find_extrema("spots/I6_PTmax.txt", "PTmax", 20.d0*coop_SI_arcmin)
  
end program test
