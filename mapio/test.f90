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

  type(coop_healpix_maps)::hgm
  call coop_healpix_spot_select_mask(nside=1024, l_deg = 207.8d0, b_deg = -56.3d0, r_deg = 12.d0, filename = "planck14/coldspot_select_mask_1024.fits")
  call coop_healpix_merge_masks("planck14/coldspot_select_mask_1024.fits", "planck14/dx11_v2_common_int_mask_010a_1024.fits", "planck14/coldspot_select_mask_1024.fits")  

end program stackth
