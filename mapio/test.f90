program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"

  character(LEN=*),parameter::mapdir = "../data/cmb/maps/act/"
  character(LEN=*),parameter::fitsfile = mapdir//"beam60.fits"
  type(coop_fits_image_cea)::cf
  type(coop_sphere_disc) disc
  COOP_INT nbins, i
  real density(1000), center(10000)
  COOP_LONG_INT pix
  COOP_INT ix, iy
  COOP_REAL coor(2)
  call cf%open(fitsfile)
  call cf%header%print()
!!$  call cf%get_data()
!!$  where(abs(cf%image) .gt. 3.e3)
!!$     cf%image = 0.
!!$  end where
!!$  call coop_bin_data(cf%image, nbins, center, density)
!!$  do i=1, nbins
!!$     write(*,*) i, center(i), density(i)
!!$  enddo
!!$  disc = coop_sphere_disc(coop_pio2, 0.d0, coop_pio2, 1.d0)
end program test
