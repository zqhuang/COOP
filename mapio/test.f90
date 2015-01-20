program stackth
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  call coop_healpix_latitude_cut_mask(256, 30.d0, "planck14/lat30_mask_n256.fits")
end program stackth
