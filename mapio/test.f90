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

  call coop_healpix_spot_select_mask(nside=1024, l_deg = 212.d0, b_deg = -13.d0, r_deg = 90.d0, filename = "planck14/power_asymmetry_hemisphere_north_1024.fits")
  call coop_healpix_spot_cut_mask(nside=1024, l_deg = 212.d0, b_deg = -13.d0, r_deg = 90.d0, filename = "planck14/power_asymmetry_hemisphere_south_1024.fits")
  
  call coop_healpix_merge_masks("planck14/power_asymmetry_hemisphere_north_1024.fits", "planck14/dx11_v2_common_int_mask_010a_1024.fits", "planck14/power_asymmetry_hemisphere_north_int_1024.fits")
  call coop_healpix_merge_masks("planck14/power_asymmetry_hemisphere_south_1024.fits", "planck14/dx11_v2_common_int_mask_010a_1024.fits", "planck14/power_asymmetry_hemisphere_south__int_1024.fits")  

  call coop_healpix_merge_masks("planck14/power_asymmetry_hemisphere_north_1024.fits", "planck14/dx11_v2_common_pol_mask_010a_1024.fits", "planck14/power_asymmetry_hemisphere_north_pol_1024.fits")
  call coop_healpix_merge_masks("planck14/power_asymmetry_hemisphere_south_1024.fits", "planck14/dx11_v2_common_pol_mask_010a_1024.fits", "planck14/power_asymmetry_hemisphere_south__pol_1024.fits")


  call coop_healpix_spot_cut_mask(nside=1024, l_deg = 207.8d0, b_deg = -56.3d0, r_deg = 6.d0, filename = "planck14/coldspot_cut_mask_1024.fits")
  
  call coop_healpix_merge_masks("planck14/coldspot_cut_mask_1024.fits", "planck14/dx11_v2_common_int_mask_010a_1024.fits", "planck14/coldspot_cut_mask_int_1024.fits")

  call coop_healpix_merge_masks("planck14/coldspot_cut_mask_1024.fits", "planck14/dx11_v2_common_pol_mask_010a_1024.fits", "planck14/coldspot_mask_pol_1024.fits")



  call coop_healpix_spot_select_mask(nside=512, l_deg = 212.d0, b_deg = -13.d0, r_deg = 90.d0, filename = "planck14/power_asymmetry_hemisphere_north_0512.fits")
  call coop_healpix_spot_cut_mask(nside=512, l_deg = 212.d0, b_deg = -13.d0, r_deg = 90.d0, filename = "planck14/power_asymmetry_hemisphere_south_0512.fits")
  
  call coop_healpix_merge_masks("planck14/power_asymmetry_hemisphere_north_0512.fits", "planck14/dx11_v2_common_int_mask_020a_0512.fits", "planck14/power_asymmetry_hemisphere_north_int_0512.fits")
  call coop_healpix_merge_masks("planck14/power_asymmetry_hemisphere_south_0512.fits", "planck14/dx11_v2_common_int_mask_020a_0512.fits", "planck14/power_asymmetry_hemisphere_south__int_0512.fits")  

  call coop_healpix_merge_masks("planck14/power_asymmetry_hemisphere_north_0512.fits", "planck14/dx11_v2_common_pol_mask_020a_0512.fits", "planck14/power_asymmetry_hemisphere_north_pol_0512.fits")
  call coop_healpix_merge_masks("planck14/power_asymmetry_hemisphere_south_0512.fits", "planck14/dx11_v2_common_pol_mask_020a_0512.fits", "planck14/power_asymmetry_hemisphere_south__pol_0512.fits")


  call coop_healpix_spot_cut_mask(nside=512, l_deg = 207.8d0, b_deg = -56.3d0, r_deg = 6.d0, filename = "planck14/coldspot_cut_mask_0512.fits")
  
  call coop_healpix_merge_masks("planck14/coldspot_cut_mask_0512.fits", "planck14/dx11_v2_common_int_mask_020a_0512.fits", "planck14/coldspot_cut_mask_int_0512.fits")

  call coop_healpix_merge_masks("planck14/coldspot_cut_mask_0512.fits", "planck14/dx11_v2_common_pol_mask_020a_0512.fits", "planck14/coldspot_cut_mask_pol_0512.fits")     
  
  
end program stackth
