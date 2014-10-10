program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"

  COOP_INT::fwhm_arcmin = 15
  COOP_STRING::spot_type = "PmaxSortT"
  COOP_STRING::input_file ="dust/dust_siqu.fits"
  COOP_STRING::imask_file  ="ffp7/ffp7_union_polmask_1024.fits"
  COOP_STRING::polmask_file  ="commander/commander_polmask.fits"
  COOP_REAL:: threshold = 1.

  COOP_STRING::mask_file
  COOP_STRING output_file
  COOP_REAL::fwhm
  COOP_INT l
  COOP_STRING prefix, force_outfile

  if(iargc().ge. 7)then
     input_file=trim(adjustl(coop_inputargs(1)))
     spot_type=trim(adjustl(coop_inputargs(2)))
     imask_file=trim(adjustl(coop_inputargs(3)))
     polmask_file=trim(adjustl(coop_inputargs(4)))
     force_outfile = trim(adjustl(coop_inputargs(5)))
     threshold = coop_str2int(coop_InputArgs(6))
     fwhm_arcmin = coop_str2int(coop_InputArgs(7))
  else
     force_outfile = ""
  endif
  fwhm = fwhm_arcmin*coop_SI_arcmin
  select case(trim(spot_type))
  case("Tmax", "Tmax_QTUTOrient", "PTmax", "Tmin", "Tmin_QTUTOrient", "PTmin", "PTmaxSortT")
     mask_file = trim(adjustl(imask_file))
  case("Pmax", "Pmin", "Emax", "Emin", "Bmax", "Bmin","PmaxSortT")
     mask_file = trim(adjustl(polmask_file))
  case default
     write(*,*) trim(spot_type)
     stop "unknown spot type"
  end select
  if(trim(force_outfile).eq."")then
     prefix = "spots/"//trim(coop_file_name_of(input_file, want_ext =  .false.))//"_"
  else
     prefix = "spots/"//trim(force_outfile)//"_"
  endif
  if(abs(threshold) .lt. 6)then
     output_file = trim(prefix)//trim(spot_type)//"_threshold"//trim(coop_num2str(nint(threshold)))//"_fwhm"//trim(coop_num2str(fwhm_arcmin))//".txt"
     call coop_healpix_export_spots(trim(input_file), trim(output_file), trim(spot_type), threshold, mask_file = trim(mask_file), filter_fwhm = fwhm)
  else
     output_file = trim(prefix)//trim(spot_type)//"_NoThreshold_fwhm"//trim(coop_num2str(fwhm_arcmin))//".txt"
     call coop_healpix_export_spots(trim(input_file), trim(output_file), trim(spot_type), mask_file = trim(mask_file), filter_fwhm = fwhm)
  endif
  write(*,*) "The output file is: "//trim(output_file)
end program test
