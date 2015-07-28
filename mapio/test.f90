program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::mask
  call mask%generate_latcut_mask(nside =1, latitude_deg = 30.d0)! #, depth_deg = 10.d0)
  call mask%write("dust/lat30_mask_n1.fits")
end program test
