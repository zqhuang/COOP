program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"

  COOP_INT::fwhm_arcmin = 15.
  COOP_INT::fwhm_in = 15.
  COOP_STRING::spot_type = "Tmax_QTUTOrient"
  COOP_STRING::input_file ="massive/simu_TQTUT_00000_015a_n1024.fits"
  COOP_STRING::imask_file  = "" 
  COOP_STRING::polmask_file  = "" 
  COOP_REAL:: threshold = 1.e31
  COOP_REAL:: threshold_pol = 1.e31

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
  fwhm = sqrt(max(dble(fwhm_arcmin)**2 - dble(fwhm_in)**2, 0.d0))*coop_SI_arcmin
  select case(trim(spot_type))
  case("Tmax", "Tmax_QTUTOrient", "PTmax", "Tmin", "Tmin_QTUTOrient", "PTmin", "PTmaxSortT", "zetamax", "zetamax_qzuzOrient", "zetamin", "zetamin_qzuzOrient", "PZmax", "PZmin")
     mask_file = trim(adjustl(imask_file))
  case("Pmax", "Pmin", "Emax", "Emin", "Bmax", "Bmin","PmaxSortT")
     mask_file = trim(adjustl(polmask_file))
  case default
     write(*,*) trim(spot_type)
     stop "unknown spot type"
  end select
  prefix = "spots/"//trim(coop_file_name_of(input_file, want_ext =  .false.))//"_fwhm"//trim(coop_num2str(fwhm_arcmin))//"_"     
  if(trim(force_outfile).ne."")then
     output_file = trim(adjustl(force_outfile))
  else
     output_file = trim(prefix)//trim(spot_type)
     if(threshold.lt.coop_healpix_max_threshold) &
          output_file = trim(output_file)//"_threshold"//trim(coop_numstr2goodstr(coop_num2str(threshold)))
     if(threshold_pol.lt.coop_healpix_max_threshold) &
          output_file = trim(output_file)//"_2ndthreshold"//trim(coop_numstr2goodstr(coop_num2str(threshold_pol)))
     output_file = trim(output_file)//".txt"
  endif
  call coop_healpix_export_spots(trim(input_file), trim(output_file), trim(spot_type), threshold = threshold, threshold_pol = threshold_pol, mask_file = trim(mask_file), fwhm = fwhm)  
  write(*,*) "The output file is: "//trim(output_file)
end program test
