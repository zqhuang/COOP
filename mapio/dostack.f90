program test
  use coop_healpix_mod
  use coop_wrapper_utils
  implicit none

#define USE_PLANCK 1

#include "constants.h"

  
  COOP_REAL, parameter::smooth_fwhm = 15.*coop_SI_arcmin
  COOP_UNKNOWN_STRING,parameter:: color_table = "Rainbow"
  COOP_UNKNOWN_STRING, parameter :: spot_type = "B"
  COOP_REAL,parameter::r=2.*coop_SI_degree, dr = max(smooth_fwhm/3., r/50.)
  COOP_INT, parameter::n = ceiling(r/dr)

#ifdef USE_PLANCK
  COOP_UNKNOWN_STRING, parameter :: map_file = "inps/predx11_iqu_nside512_inp_teb0200_submap003.fits"
  COOP_UNKNOWN_STRING, parameter :: spots_file = "spots/predx11_iqu_nside512_inp_teb0200_submap003_Bmax_threshold0_fwhm15.txt" 
  COOP_UNKNOWN_STRING, parameter :: imask_file = "inps/predx11_imask_nside512.fits" 
  COOP_UNKNOWN_STRING, parameter :: polmask_file = "inps/predx11_polmask_nside512.fits"
#endif

#ifdef USE_WMAP
  COOP_UNKNOWN_STRING, parameter :: map_file = "wmap/wmap_iqu.fits"
  COOP_UNKNOWN_STRING, parameter :: spots_file = "spots/wmap_iqu_Tmax_NoThreshold_fwhm30.txt"
  COOP_UNKNOWN_STRING, parameter :: imask_file = "wmap/wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits" 
  COOP_UNKNOWN_STRING, parameter :: polmask_file = "wmap/wmap_polarization_analysis_mask_r9_9yr_v5.fits"
#endif

#ifdef USE_SIMU
  COOP_UNKNOWN_STRING, parameter :: map_file = "simu/sim3_T.fits" 
  COOP_UNKNOWN_STRING, parameter :: spots_file = "spots/sim3_T_Tmax_QTUTOrient_NoThreshold_fwhm10.txt"
  COOP_UNKNOWN_STRING, parameter :: imask_file = ""
  COOP_UNKNOWN_STRING, parameter :: polmask_file = ""
#endif

#ifdef USE_FFP7
  COOP_UNKNOWN_STRING, parameter :: map_file = "ffp7/ffp7_iqu_smoothed_fwhm15arcmin.fits" 
  COOP_UNKNOWN_STRING, parameter :: spots_file = "spots/ffp7_TQTUT_TQUmax_NoThreshold_fwhm15.txt" 
  COOP_UNKNOWN_STRING, parameter :: imask_file = "ffp7/ffp7_smica_harmonic_mask_I.fits"
  COOP_UNKNOWN_STRING, parameter :: polmask_file = "ffp7/ffp7_smica_harmonic_mask_QU.fits"
#endif

#ifdef USE_PREDX11
  COOP_UNKNOWN_STRING, parameter :: map_file = "predx11/predx11_iqu_smoothed_fwhm15arcmin.fits" 
  COOP_UNKNOWN_STRING, parameter :: spots_file = "spots/predx11_TQTUT_TQUmax_NoThreshold_fwhm15.txt" 
  COOP_UNKNOWN_STRING, parameter :: imask_file = "predx11/predx11_imask.fits"
  COOP_UNKNOWN_STRING, parameter :: polmask_file = "predx11/predx11_polmask.fits"
#endif

  COOP_UNKNOWN_STRING, parameter :: prefix = "stacked/"
  COOP_STRING fout, caption, fname
  COOP_INT,parameter::mmax = 4
  type(coop_healpix_maps) map, mask
  type(coop_healpix_patch) patch
  logical::do_mask
  select case(trim(spot_type))
  case("T", "I")
     if(trim(imask_file) .ne. "")then
        do_mask = .true.
        call mask%read( imask_file, nmaps_wanted = 1 )
     else
        do_mask = .false.
     endif
  case("E","B","Q", "QrUr","QU")
     if(trim(polmask_file) .ne. "")then
        call mask%read( polmask_file, nmaps_wanted = 1)
        do_mask = .true.
     else
        do_mask = .false.
     endif
  case default
     stop "Unknown spot_type"
  end select
  call map%read(map_file)

  if(index(spots_file, "_Tmax_QTUTOrient_") .gt. 0 .or. index(spots_file, "TQUmax").gt. 0)then
     caption = "$T$ maxima, oriented"
  elseif(index(spots_file, "_Tmax_").gt.0)then
     caption = "$T$ maxima, random orientation"
  elseif(index(spots_file, "_Emax_").gt.0)then
     caption = "$E$ maxima, random orientation"
  elseif(index(spots_file, "_Bmax_").gt.0)then
     caption = "$B$ maxima, random orientation"
  elseif(index(spots_file, "_PTmax_").gt.0)then
     caption = "$P_T$ maxima, oriented"
  elseif(index(spots_file, "_Pmax_").gt.0)then
     caption = "$P$ maxima, oriented"
  elseif(index(spots_file, "_Tmin_QTUTOrient_").gt.0)then
     caption = "$T$ minima, oriented"
  elseif(index(spots_file, "_Tmin_") .gt. 0)then
     caption = "$T$ minima, random orientation"
  elseif(index(spots_file, "_PTmin_").gt.0)then
     caption = "$P_T$ minima, oriented"
  elseif(index(spots_file, "_Pmin_").gt.0)then
     caption = "$P$ minima, oriented"
  else
     stop "Unknown spots_file class"
  endif

  if(index(spots_file, "threshold0").gt.0)then
     caption = trim(caption)//", threshold $\nu$=0"
  elseif(index(spots_file, "threshold1").gt.0)then
     caption = trim(caption)//", threshold $\nu$=1"
  elseif(index(spots_file, "threshold2").gt.0)then
     caption = trim(caption)//", threshold $\nu$=2"
  elseif(index(spots_file, "threshold3").gt.0)then
     caption = trim(caption)//", threshold $\nu$=3"
  elseif(index(spots_file, "threshold4").gt.0)then
     caption = trim(caption)//", threshold $\nu$=4"
  elseif(index(spots_file, "threshold5").gt.0)then
     caption = trim(caption)//", threshold $\nu$=5"
  endif

  call patch%init(spot_type, n, dr, mmax = mmax)
  if(do_mask)then
     call map%stack(patch, spots_file, mask)
  else
     call map%stack(patch, spots_file)
  endif
  caption = trim(coop_num2str(patch%nstack_raw))//" patches on "//trim(caption)
  fname = coop_file_name_of(spots_file)
  select case(spot_type)
  case("QrUr")
     fout = prefix//"Qr_on_"//trim(fname)
     call patch%plot(imap = 1, output =trim(fout), label = "$Q_r(\mu K)$", caption = caption, color_table = color_table, headless_vectors = .true.)
     call patch%plot(imap = 2, output = prefix//"Ur_on_"//trim(fname), label = "$U_r(\mu K)$", caption = caption, color_table = color_table, headless_vectors = .true.)
  case("QU")
     fout = prefix//"Q_on_"//trim(fname)
     call patch%plot(imap = 1, output = trim(fout), label = "$Q(\mu K)$", caption = caption, color_table = color_table, headless_vectors = .true.)
     call patch%plot(imap = 2, output = prefix//"U_on_"//trim(fname), label = "$U(\mu K)$", caption = caption, color_table = color_table, headless_vectors = .true.)
  case("T", "E", "B", "I") 
     fout = prefix//spot_type//"_on_"//trim(fname)
     call patch%plot(imap = 1, output =trim(fout), label = "$"//spot_type//"(\mu K)$", caption = caption, color_table = color_table, headless_vectors = .true.)
  end select
  write(*,*) "the output file is: "//trim(fout)
end program test
