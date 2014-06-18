program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"

  integer,parameter::fwhm_arcmin = 15
  COOP_UNKNOWN_STRING, parameter::spot_type = "PTmax"
  COOP_UNKNOWN_STRING, parameter::input_file ="planck/smica_inp_cmb.fits" 
  COOP_UNKNOWN_STRING, parameter::mask_file  = "planck/smica_valmask.fits"  

  COOP_REAL, parameter:: threshold = 1.e31
  COOP_STRING output_file
  COOP_REAL,parameter::fwhm = fwhm_arcmin*coop_SI_arcmin
  COOP_INT l
  COOP_STRING prefix

  prefix = "spots/"//trim(coop_file_name_of(input_file, want_ext =  .false.))//"_"
  print*, " filter fwhm = ", fwhm_arcmin , " arcmin"
  if(abs(threshold).lt. 1.e20)then
     output_file = trim(prefix)//trim(spot_type)//"_threshold"//trim(coop_num2str(threshold))//"_fwhm"//trim(coop_num2str(fwhm_arcmin))//".txt"
     call coop_healpix_export_spots(input_file, trim(output_file), spot_type, threshold, mask_file = mask_file, filter_fwhm = fwhm)
  else
     output_file = trim(prefix)//trim(spot_type)//"_NoThreshold_fwhm"//trim(coop_num2str(fwhm_arcmin))//".txt"
     call coop_healpix_export_spots(input_file, trim(output_file), spot_type, mask_file = mask_file, filter_fwhm = fwhm)
  endif
end program test
