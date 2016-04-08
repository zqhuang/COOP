program shells
  use coop_wrapper_firstorder
  use coop_zeta3d_mod
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::mask
  call mask%generate_latcut_mask(nside = 512, latitude_deg = 30.d0, depth_deg = 1.d0)
  call mask%write("planck15/mask_lat30_n512.fits")
end program shells
