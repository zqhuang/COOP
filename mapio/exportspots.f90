program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"

  COOP_INT::fwhm_arcmin = 0.
  COOP_INT::fwhm_in = 0.
  COOP_STRING::spot_type = "Bmax"
  COOP_STRING::input_file ="simu/simurp2_iqu_10arc_n1024_converted_to_TEB_submap003.fits"
  COOP_STRING::imask_file  ="commander/commander_imask.fits"
  COOP_STRING::polmask_file  ="commander/commander_polmask.fits"
  COOP_REAL:: threshold = 0.

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
     fwhm_in = coop_str2int(coop_InputArgs(8))
  else
     force_outfile = ""
  endif
  fwhm = sqrt(dble(fwhm_arcmin)**2 - dble(fwhm_in)**2)*coop_SI_arcmin
  select case(trim(spot_type))
  case("Tmax", "Tmax_QTUTOrient", "PTmax", "Tmin", "Tmin_QTUTOrient", "PTmin", "PTmaxSortT", "zetamax")
     mask_file = trim(adjustl(imask_file))
  case("Pmax", "Pmin", "Emax", "Emin", "Bmax", "Bmin","PmaxSortT")
     mask_file = trim(adjustl(polmask_file))
  case default
     write(*,*) trim(spot_type)
     stop "unknown spot type"
  end select
  prefix = "spots/"//trim(coop_file_name_of(input_file, want_ext =  .false.))//"_fwhm"//trim(coop_num2str(fwhm_arcmin))//"_"     
  if(abs(threshold) .lt. 6)then
     if(trim(force_outfile).ne."")then
        output_file = trim(adjustl(force_outfile))
     else
        output_file = trim(prefix)//trim(spot_type)//"_threshold"//trim(coop_num2str(nint(threshold)))//".txt"
     endif
     call coop_healpix_export_spots(trim(input_file), trim(output_file), trim(spot_type), threshold = threshold, mask_file = trim(mask_file), fwhm = fwhm)
  else
     if(trim(force_outfile).ne."")then
        output_file = trim(adjustl(force_outfile))
     else
        output_file = trim(prefix)//trim(spot_type)//"_NoThreshold.txt"
     endif
     call coop_healpix_export_spots(trim(input_file), trim(output_file), trim(spot_type), mask_file = trim(mask_file), fwhm = fwhm)
  endif
  write(*,*) "The output file is: "//trim(output_file)
end program test
