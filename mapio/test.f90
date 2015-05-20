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
  type(coop_file)::fp
  COOP_STRING::line
  COOP_INT ::i
  call map%read("planck14/dx11_v2_smica_pol_case1_cmb_hp_20_40_010a_1024.fits")
  call mask%read("planck14/dx11_v2_smica_pol_mask_010a_1024.fits")
  map%map = map%map*1.e6
  where(mask%map(:,1).eq.0.)
     map%map(:,1) = -1.63750e30
     map%map(:,2) = -1.63750e30          
  end where
  call map%write("planck14/smica_masked.fits")

  
end program test
