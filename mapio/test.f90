program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::mask
  COOP_INT::l
  call mask%init(nside = 128, nmaps = 1, genre  = "MASK")
  call coop_healpix_mask_hemisphere(mask, 0.d0, 0.d0)
  call mask%write("lowl/hemisphere_mask.fits")
end program test
