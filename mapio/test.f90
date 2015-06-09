program test
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
!!$  call coop_healpix_latitude_cut_mask(256, 30.d0, "lowl/lat30_mask_n0256.fits")  
!!$  call coop_healpix_latitude_cut_mask(256, 20.d0, "lowl/lat20_mask_n0256.fits")
!!$  call coop_healpix_latitude_cut_mask(256, 10.d0, "lowl/lat10_mask_n0256.fits")
!!$  call coop_healpix_latitude_cut_mask(256, 15.d0, "lowl/lat15_mask_n0256.fits")  
!!$  call coop_healpix_latitude_cut_mask(256, 5.d0, "lowl/lat5_mask_n0256.fits")
!!$  call coop_healpix_latitude_cut_mask(256, 3.d0, "lowl/lat3_mask_n0256.fits")
!!$  stop
!!$  call coop_healpix_mask_reverse("lowl/lat30_mask_n0256.fits", "lowl/lat30_revmask_n0256.fits")  
!!$  call coop_healpix_mask_reverse("lowl/lat20_mask_n0256.fits", "lowl/lat20_revmask_n0256.fits")
 call coop_healpix_mask_reverse("lowl/lat15_mask_n0256.fits", "lowl/lat15_revmask_n0256.fits")  
 !!$ call coop_healpix_mask_reverse("lowl/lat10_mask_n0256.fits", "lowl/lat10_revmask_n0256.fits")
!!$ call coop_healpix_mask_reverse("lowl/lat5_mask_n0256.fits", "lowl/lat5_revmask_n0256.fits")
!!$  call coop_healpix_mask_reverse("lowl/lat3_mask_n0256.fits", "lowl/lat3_revmask_n0256.fits")    
end program test
