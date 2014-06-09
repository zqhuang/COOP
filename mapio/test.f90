program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"

  character(LEN=*),parameter::mapdir = "../data/cmb/maps/act/"
  character(LEN=*),parameter::fitsfile = mapdir//"Q50.fits"
  type(coop_fits_image_cea)::cf
  type(coop_sphere_disc) disc
  COOP_LONG_INT pix
  COOP_INT ix, iy
  COOP_REAL coor(2)
  call cf%open(fitsfile)
  print*, "nside = ", cf%nside
  print*, "npix = ", cf%npix
  call cf%get_data()
  disc = coop_sphere_disc(coop_pio2, 0.d0, coop_pio2, 1.d0)
  do ix = -2, 2
     do iy = -2, 2
        pix = (cf%center(2)-1 + iy)*cf%nside(1) + cf%center(1) - 1+ix
        call cf%pix2flat(disc, pix,  coor)
        write(*,*) ix, iy, coor
     enddo
  enddo
end program test
