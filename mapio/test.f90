program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  call coop_healpix_latitude_cut_smoothmask(nside = 1024, latitude_degree = 30.d0, depth = 5.d0, filename = "tuhin/lat30_mask_n1024_smooth.fits")
end program test
