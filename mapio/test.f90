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
  type(coop_healpix_maps)::hgm, mask
  call hgm%read("planck14/", nmaps_wanted = 3, nmaps_to_read = 1, spin=(/ 0, 2, 2 /) )
  call hgm%iqu2TQTUT1()
  call hgm%write("planck14/dx11_v2_smica_TQTUT_cmb_010a_1024_hp_10_20.fits")
end program stackth
