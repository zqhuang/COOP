program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  use udgrade_nr
  use coord_v_convert,only:coordsys2euler_zyz
#endif  
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask
  call map%read("tuhin/dust_TQUL_015a_b30-500_n512.fits", nmaps_wanted = 4)
  call mask%read("planck14/lat30_mask_n512.fits")
  call map%distr_nu_e(mask = mask, numin = -5.d0, numax = 8.d0, emin = 0.d0, emax = 2.d0, nnu = 50, ne = 50, output = "CMB_nue_distr.txt")

  
end program test
