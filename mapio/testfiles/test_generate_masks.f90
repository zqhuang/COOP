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
  COOP_INT::nside, fwhm
  COOP_STRING::postfix, hastr, cstr
  nside = 256
  hastr = "HemAsym_"
  cstr = "ColdSpotCut_"
  do while (nside .lt. 1024)
     nside = nside*2
     fwhm = 10240/nside
     postfix = "_"//trim(coop_ndigits(fwhm, 3))//"a_"//trim(coop_ndigits(nside,4))//".fits"
     
     call coop_healpix_spot_select_mask(nside=nside, l_deg = 212.d0, b_deg = -13.d0, r_deg = 90.d0, filename = "planck14/"//trim(hastr)//"north_mask"//trim(postfix))
     call coop_healpix_spot_cut_mask(nside=nside, l_deg = 212.d0, b_deg = -13.d0, r_deg = 90.d0, filename = "planck14/"//trim(hastr)//"south_mask"//trim(postfix))

     call coop_healpix_merge_masks("planck14/"//trim(hastr)//"north_mask"//trim(postfix), "planck14/dx11_v2_common_int_mask"//trim(postfix), "planck14/"//trim(hastr)//"north_int_mask"//trim(postfix))
     call coop_healpix_merge_masks("planck14/"//trim(hastr)//"north_mask"//trim(postfix), "planck14/dx11_v2_common_pol_mask"//trim(postfix), "planck14/"//trim(hastr)//"north_pol_mask"//trim(postfix))     

     call coop_healpix_merge_masks("planck14/"//trim(hastr)//"south_mask"//trim(postfix), "planck14/dx11_v2_common_int_mask"//trim(postfix), "planck14/"//trim(hastr)//"south_int_mask"//trim(postfix))
     call coop_healpix_merge_masks("planck14/"//trim(hastr)//"south_mask"//trim(postfix), "planck14/dx11_v2_common_pol_mask"//trim(postfix), "planck14/"//trim(hastr)//"south_pol_mask"//trim(postfix))     

     


     call coop_healpix_spot_cut_mask(nside=nside, l_deg = 207.8d0, b_deg = -56.3d0, r_deg = 6.d0, filename = "planck14/"//trim(cstr)//"mask"//trim(postfix))

     call coop_healpix_merge_masks("planck14/"//trim(cstr)//"mask"//trim(postfix), "planck14/dx11_v2_common_int_mask"//trim(postfix), "planck14/"//trim(cstr)//"int_mask"//trim(postfix) )

     call coop_healpix_merge_masks("planck14/"//trim(cstr)//"mask"//trim(postfix), "planck14/dx11_v2_common_pol_mask"//trim(postfix), "planck14/"//trim(cstr)//"pol_mask"//trim(postfix) )     

  enddo


  
end program stackth
