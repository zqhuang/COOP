program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::hm
  type(coop_fits_image_cea)::fm
  call fm%read("act16/deep56_coadd_I.fits")
  call hm%read("planck14/COM_CMB_IQU-smica_1024_R2.02_full.fits", nmaps_wanted = 3)
  hm%map = hm%map*1.d6
  call hm%rotate_coor(from = "G", to = "C")
  call fm%from_healpix(hm, 2)
  call fm%write("act16/planck_smica_Q.fits")
  call fm%from_healpix(hm, 3)
  call fm%write("act16/planck_smica_U.fits")
end program test
