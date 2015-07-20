program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  call coop_healpix_latitude_cut_smoothmask(nside = 1024, latitude_degree = 25.d0,  depth = 10.d0, filename = "dust/lat25_mask_n1024_smooth.fits")
end program test
