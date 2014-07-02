program test
  use coop_healpix_mod
  use coop_wrapper_utils
  implicit none

#define USE_PLANCK 1

#include "constants.h"

  
  COOP_REAL, parameter::pre_smooth_fwhm = 15.*coop_SI_arcmin
  COOP_UNKNOWN_STRING,parameter:: color_table = "Planck"
  COOP_UNKNOWN_STRING, parameter :: spot_type = "Qr"
  COOP_REAL ,parameter::r=2.*coop_SI_degree, dr = max(coop_SI_arcmin/2., r/30.)

#ifdef USE_PLANCK
  COOP_UNKNOWN_STRING, parameter :: map_file = "predx11/predx11_iqu.fits"
  COOP_UNKNOWN_STRING, parameter :: spot_file = "spots/predx11_iqu_Tmax_threshold0_fwhm15.txt" 
  COOP_UNKNOWN_STRING, parameter :: imask_file = "predx11/predx11_sk.fits" 
  COOP_UNKNOWN_STRING, parameter :: polmask_file = "predx11/predx11_polmask.fits"
#endif

#ifdef USE_WMAP
  COOP_UNKNOWN_STRING, parameter :: map_file = "wmap/wmap_iqu.fits"
  COOP_UNKNOWN_STRING, parameter :: spot_file = "spots/wmap_iqu_Tmax_NoThreshold_fwhm30.txt"
  COOP_UNKNOWN_STRING, parameter :: imask_file = "wmap/wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits" 
  COOP_UNKNOWN_STRING, parameter :: polmask_file = "wmap/wmap_polarization_analysis_mask_r9_9yr_v5.fits"
#endif

#ifdef USE_SIMU
  COOP_UNKNOWN_STRING, parameter :: map_file = "simu/sim3_T.fits" 
  COOP_UNKNOWN_STRING, parameter :: spot_file = "spots/sim3_T_Tmax_QTUTOrient_NoThreshold_fwhm10.txt"
  COOP_UNKNOWN_STRING, parameter :: imask_file = ""
  COOP_UNKNOWN_STRING, parameter :: polmask_file = ""
#endif

#ifdef USE_FFP7
  COOP_UNKNOWN_STRING, parameter :: map_file = "ffp7/ffp7_iqu_smoothed_fwhm15arcmin.fits" 
  COOP_UNKNOWN_STRING, parameter :: spot_file = "spots/ffp7_TQTUT_TQUmax_NoThreshold_fwhm15.txt" 
  COOP_UNKNOWN_STRING, parameter :: imask_file = "ffp7/ffp7_smica_harmonic_mask_I.fits"
  COOP_UNKNOWN_STRING, parameter :: polmask_file = "ffp7/ffp7_smica_harmonic_mask_QU.fits"
#endif

#ifdef USE_PREDX11
  COOP_UNKNOWN_STRING, parameter :: map_file = "predx11/predx11_iqu_smoothed_fwhm15arcmin.fits" 
  COOP_UNKNOWN_STRING, parameter :: spot_file = "spots/predx11_TQTUT_TQUmax_NoThreshold_fwhm15.txt" 
  COOP_UNKNOWN_STRING, parameter :: imask_file = "predx11/predx11_imask.fits"
  COOP_UNKNOWN_STRING, parameter :: polmask_file = "predx11/predx11_polmask.fits"
#endif

  COOP_UNKNOWN_STRING, parameter :: output_prefix = "stacked/stack_"
  COOP_STRING spotname
  COOP_STRING mask_file
  COOP_STRING fout, caption

  select case(trim(spot_type))
  case("T")
     mask_file = imask_file
  case("E","B","Q", "Qr","U", "Ur")
     mask_file = polmask_file
  case default
     mask_file = ""
     stop "Unknown spot_type"
  end select
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
        call coop_healpix_stack_io(map_file, trim(fout),  spot_file, r, dr, "Qr", "$Q^T_{r}(\mu K)$", headless_vector=.true., caption = trim(caption), mask_file = trim(mask_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table)
     else
        fout = trim(output_prefix)//"Qr_"//trim(spotname)
        call coop_healpix_stack_io(map_file, trim(fout), spot_file, r, dr, "Qr", "$Q_r(\mu K)$",  headless_vector=.true., caption = trim(caption), mask_file = trim(mask_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table ) 
     endif
  case("Ur")
     if(index(map_file, "_QTUT").gt.0)then
        fout = trim(output_prefix)//"UTr_"//trim(spotname)
        call coop_healpix_stack_io(map_file, trim(fout),  spot_file, r, dr, "Ur", "$U^T_r(\mu K)$", headless_vector=.true., caption = trim(caption), mask_file = trim(mask_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table ) 

     else
        fout = trim(output_prefix)//"Ur_"//trim(spotname)
        call coop_healpix_stack_io(map_file, trim(fout),  spot_file, r, dr, "Ur", "$U_r(\mu K)$", headless_vector=.true., caption = trim(caption), mask_file = trim(mask_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table ) 
     endif
  case("Q") 
     if(index(map_file, "_QTUT").gt.0)then
        fout = trim(output_prefix)//"QT_"//trim(spotname)
        call coop_healpix_stack_io(map_file, trim(fout), spot_file, r, dr, "Q", "$Q^T(\mu K)$", headless_vector=.true. , m_filter=(/ 0, 2, 4 /), caption = trim(caption), mask_file = trim(mask_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table )
     else    
        fout = trim(output_prefix)//"Q_"//trim(spotname)
        call coop_healpix_stack_io(map_file, trim(fout), spot_file, r, dr, "Q", "$Q(\mu K)$", headless_vector=.true. , m_filter=(/ 0, 2, 4 /), caption = trim(caption), mask_file = trim(mask_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table )
     endif
  case("U")
     if(index(map_file, "_QTUT").gt.0)then
        fout = trim(output_prefix)//"UT_"//trim(spotname)
        call coop_healpix_stack_io(map_file, trim(fout), spot_file, r, dr, trim(spot_type), "$U^T(\mu K)$", headless_vector=.true. , caption = trim(caption), mask_file = trim(mask_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table)

     else
        fout = trim(output_prefix)//"U_"//trim(spotname)
        call coop_healpix_stack_io(map_file, trim(fout), spot_file, r, dr, trim(spot_type), "$U(\mu K)$", headless_vector=.true. , caption = trim(caption), mask_file = trim(mask_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table)
     endif
  case("T","E","B")
     fout = trim(output_prefix)//trim(spot_type)//"_"//trim(spotname)
     call coop_healpix_stack_io(map_file, trim(fout), spot_file, r, dr, trim(spot_type), "$"//trim(spot_type)//"(\mu K)$", headless_vector=.true. , caption = trim(caption), mask_file = trim(mask_file), pre_smooth_fwhm = pre_smooth_fwhm, color_table = color_table)
  case default
     stop "Unknown spot_type"
  end select
  write(*,*) "the output file is: "//trim(fout)
end program test
