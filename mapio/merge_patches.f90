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
  type(coop_healpix_patch)::hot, cold
  
  call hot%import("stacked/commander_mask94_hotT_Randrot_n0016_440a.patch")
  call cold%import("stacked/commander_mask94_coldT_Randrot_n0016_440a.patch")
  hot%image(:,:,1) = hot%image(:,:,1)*hot%nstack - cold%image(:,:,1)*cold%nstack
  hot%nstack = hot%nstack + cold%nstack
  hot%caption = "hot - cold"
  hot%tbs%zmax(1) = 20.
  hot%tbs%zmin(1) = -20.
  where(hot%nstack .ne. 0.d0)
     hot%image(:,:,1) = hot%image(:,:,1)/hot%nstack
  elsewhere
     hot%image(:,:,1) = 0.d0
  end where
  call hot%plot(1, "stacked/commander_mask94_hot_minus_cold_Randrot_n0016_440a.txt")
end program test
