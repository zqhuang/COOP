program test
  use coop_healpix_mod
  use coop_wrapper_utils
  implicit none


#include "constants.h"

  
  COOP_REAL, parameter::smooth_fwhm = 15.*coop_SI_arcmin
  COOP_UNKNOWN_STRING,parameter:: color_table = "Planck"
  COOP_UNKNOWN_STRING, parameter :: spot_type = "QU"
  COOP_REAL,parameter::r=7.*coop_SI_degree, dr = max(smooth_fwhm/2., r/100.)
  COOP_INT, parameter::n = ceiling(r/dr)

  COOP_INT, parameter::npatches = 10.
  real, parameter::uppercut  = 115.
  real, parameter::lowercut  = 15.  
  COOP_REAL,parameter::zmin = -1.
  COOP_REAL,parameter::zmax = 8.
!!$  COOP_UNKNOWN_STRING, parameter :: map_file = "planck14/commander_siqu.fits"
!!$  COOP_UNKNOWN_STRING, parameter :: spots_file = "spots/commander_siqu_fwhm15_PmaxSortT_threshold1.txt"
!!$  COOP_UNKNOWN_STRING, parameter :: prefix = "multstacked/cmb"
  
  COOP_UNKNOWN_STRING, parameter :: map_file = "dust/dust_siqu.fits"
  COOP_UNKNOWN_STRING, parameter :: spots_file = "spots/dust_siqu_fwhm15_PmaxSortT_threshold1.txt"
  COOP_UNKNOWN_STRING, parameter :: prefix = "multstacked/dust"
  
  COOP_UNKNOWN_STRING, parameter :: imask_file = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_UNKNOWN_STRING, parameter :: polmask_file = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"

  COOP_STRING fout, fname, caption_raw
  COOP_INT,parameter::mmax = 4
  type(coop_healpix_maps) map, mask
  type(coop_healpix_patch) patch(npatches)
  logical::do_mask
  COOP_INT i

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
  call map%smooth(smooth_fwhm)
  if(index(spots_file, "_Tmax_QTUTOrient_") .gt. 0 .or. index(spots_file, "TQUmax").gt. 0)then
     caption_raw = "$T$ maxima, oriented"
  elseif(index(spots_file, "_Tmax_").gt.0)then
     caption_raw = "$T$ maxima, random orientation"
  elseif(index(spots_file, "_Emax_").gt.0)then
     caption_raw = "$E$ maxima, random orientation"
  elseif(index(spots_file, "_Bmax_").gt.0)then
     caption_raw = "$B$ maxima, random orientation"
  elseif(index(spots_file, "_PTmax_").gt.0)then
     caption_raw = "$P_T$ maxima, oriented"
  elseif(index(spots_file, "_Pmax").gt.0)then
     caption_raw = "$P$ maxima, oriented"
  elseif(index(spots_file, "_Tmin_QTUTOrient_").gt.0)then
     caption_raw = "$T$ minima, oriented"
  elseif(index(spots_file, "_Tmin_") .gt. 0)then
     caption_raw = "$T$ minima, random orientation"
  elseif(index(spots_file, "_PTmin_").gt.0)then
     caption_raw = "$P_T$ minima, oriented"
  elseif(index(spots_file, "_Pmin_").gt.0)then
     caption_raw = "$P$ minima, oriented"
  else
     stop "Unknown spots_file class"
  endif

  if(index(spots_file, "threshold0.5").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu=0.5$"
  elseif(index(spots_file, "threshold1.5").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu=1.5$"
  elseif(index(spots_file, "threshold2.5").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu=2.5$"
  elseif(index(spots_file, "threshold0").gt.0)then
     if(index(spots_file, "_Pmax") .ne. 0 .and. index(spots_file, "_PTmax") .ne. 0)then
        caption_raw = trim(caption_raw)//", threshold $\nu=0$"
     endif
  elseif(index(spots_file, "threshold1").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu=1$"
  elseif(index(spots_file, "threshold2").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu=2$"
  elseif(index(spots_file, "threshold3").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu=3$"
  endif

  call patch(1)%init(spot_type, n, dr, mmax = mmax)
  patch(1)%color_table = trim(color_table)
  patch(1)%zmin = zmin
  patch(1)%zmax = zmax
  
  do i=2, npatches
     patch(i) = patch(1)
  enddo
  if(do_mask)then
     call map%multstack(npatches, patch, lowercut, uppercut, spots_file, mask)
  else
     call map%multstack(npatches, patch, lowercut, uppercut, spots_file)
  endif



  do i=1, npatches
     if(i.gt.1)then
        patch(i)%image = (patch(i-1)%image*patch(i-1)%nstack_raw + patch(i)%image*patch(i)%nstack_raw)/(patch(i-1)%nstack_raw + patch(i)%nstack_raw)
        patch(i)%nstack_raw = patch(i-1)%nstack_raw + patch(i)%nstack_raw
     endif
     patch(i)%caption = trim(coop_num2str(patch(i)%nstack_raw))//" patches on "//trim(caption_raw)//", "//trim(patch(i)%caption)
     
     fname = coop_file_name_of(spots_file)
     select case(spot_type)
     case("QrUr")
        fout = prefix//trim(coop_ndigits(i,2))//"_Qr_on_"//trim(fname)
        call patch(i)%plot(imap = 1, output =trim(fout))
        call patch(i)%plot(imap = 2, output = prefix//trim(coop_num2str(i))//"_Ur_on_"//trim(fname))
     case("QU")
        fout = prefix//trim(coop_ndigits(i,2))//"_Q_on_"//trim(fname)
        call patch(i)%plot(imap = 1, output = trim(fout))
        call patch(i)%plot(imap = 2, output = prefix//trim(coop_num2str(i))//"_U_on_"//trim(fname))
     case("T", "E", "B", "I") 
        fout = prefix//trim(coop_ndigits(i,2))//"_"//spot_type//"_on_"//trim(fname)
        call patch(i)%plot(imap = 1, output =trim(fout))
     end select
     write(*,*) "the output file is: "//trim(fout)
  enddo
end program test
