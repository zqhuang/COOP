program test
  use coop_healpix_mod
  use coop_wrapper_utils
  implicit none

#include "constants.h"

  
  COOP_REAL, parameter :: pre_smooth_fwhm = 15.*coop_SI_arcmin
  COOP_UNKNOWN_STRING, parameter :: color_table = "Rainbow"
  COOP_UNKNOWN_STRING, parameter :: spot_type = "T"
  COOP_REAL, parameter::r=2.*coop_SI_degree, dr = max(coop_SI_arcmin/2., r/40.)
  COOP_UNKNOWN_STRING, parameter :: map_file = "predx11/predx11_iqu.fits"
  COOP_UNKNOWN_STRING, parameter :: spot_file = "spots/predx11_iqu_Tmax_threshold0_fwhm15.txt"
!  COOP_UNKNOWN_STRING, parameter :: map_file = "ffp7/ffp7_smica_harmonic_cmb_I.fits"
! COOP_UNKNOWN_STRING, parameter :: spot_file = "spots/ffp7_smica_harmonic_cmb_I_Tmax_threshold0_fwhm15.txt"

  COOP_UNKNOWN_STRING, parameter :: mask_north_file =  "predx11/predx11_imask_north.fits"
  COOP_UNKNOWN_STRING, parameter :: mask_south_file =  "predx11/predx11_imask_south.fits"
  COOP_UNKNOWN_STRING, parameter :: output_prefix = "hsout/stack_"

  COOP_STRING spotname, label, fout, caption

  spotname = coop_file_name_of(spot_file)
  if(index(spot_file, "_Tmax_QTUTOrient_") .gt. 0 .or. index(spot_file, "TQUmax").gt. 0)then
     caption = "#$T$ maxima, oriented"
  elseif(index(spot_file, "_Tmax_").gt.0)then
     caption = "#$T$ maxima, random orientation"
  elseif(index(spot_file, "_PTmax_").gt.0)then
     caption = "#$P_T$ maxima, oriented"
  elseif(index(spot_file, "_Pmax_").gt.0)then
     caption = "#$P$ maxima, oriented"
  elseif(index(spot_file, "_Tmin_QTUTOrient_").gt.0)then
     caption = "#$T$ minima, oriented"
  elseif(index(spot_file, "_Tmin_") .gt. 0)then
     caption = "#$T$ minima, random orientation"
  elseif(index(spot_file, "_PTmin_").gt.0)then
     caption = "#$P_T$ minima, oriented"
  elseif(index(spot_file, "_Pmin_").gt.0)then
     caption = "#$P$ minima, oriented"
  else
     stop "Unknown spot_file class"
  endif
  if(index(spot_file, "threshold0").gt.0)then
     caption = trim(caption)//", threshold $\nu$=0"
  elseif(index(spot_file, "threshold1").gt.0)then
     caption = trim(caption)//", threshold $\nu$=1"
  elseif(index(spot_file, "threshold2").gt.0)then
     caption = trim(caption)//", threshold $\nu$=2"
  elseif(index(spot_file, "threshold3").gt.0)then
     caption = trim(caption)//", threshold $\nu$=3"
  elseif(index(spot_file, "threshold4").gt.0)then
     caption = trim(caption)//", threshold $\nu$=4"
  elseif(index(spot_file, "threshold5").gt.0)then
     caption = trim(caption)//", threshold $\nu$=5"
  endif
  select case(spot_type)
  case("Qr")
     if(index(map_file, "_QTUT").gt.0)then
        fout = trim(output_prefix)//"QTr_"//trim(spotname)
        label =  "$Q^T_{r}(\mu K)$"
     else
        fout = trim(output_prefix)//"Qr_"//trim(spotname)
        label =  "$Q_r(\mu K)$"
     endif
  case("Ur")
     if(index(map_file, "_QTUT").gt.0)then
        fout = trim(output_prefix)//"UTr_"//trim(spotname)
        label =  "$U_{T, r}(\mu K)$"
     else
        fout = trim(output_prefix)//"Ur_"//trim(spotname)
        label = "$U_r(\mu K)$"
     endif
  case("Q") 
     if(index(map_file, "_QTUT").gt.0)then
        fout = trim(output_prefix)//"QT_"//trim(spotname)
        label = "$Q_T(\mu K)$"
     else    
        fout = trim(output_prefix)//"Q_"//trim(spotname)
        label = "$Q(\mu K)$"
     endif
  case("U")
     if(index(map_file, "_QTUT").gt.0)then
        fout = trim(output_prefix)//"UT_"//trim(spotname)
        label  =  "$U_T(\mu K)$"

     else
        fout = trim(output_prefix)//"U_"//trim(spotname)
        label = "$U(\mu K)$"
     endif
  case("T","E","B")
     fout = trim(output_prefix)//trim(spot_type)//"_"//trim(spotname)
     label =  "$"//trim(spot_type)//"(\mu K)$"
  case default
     stop "Unknown spot_type"
  end select
  call coop_healpix_stack_io(map_file, trim(coop_file_add_postfix(fout,"_north")), spot_file, r, dr, trim(spot_type),label, headless_vector=.true. , caption = trim(caption), mask_file = trim(mask_north_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table)
  call coop_healpix_stack_io(map_file, trim(coop_file_add_postfix(fout,"_south")), spot_file, r, dr, trim(spot_type),label, headless_vector=.true. , caption = trim(caption), mask_file = trim(mask_south_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table)
  write(*,*) "the output file is: "//trim(fout)
end program test
