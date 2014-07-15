program test
  use coop_healpix_mod
  use coop_wrapper_utils
  implicit none


#include "constants.h"

  
  COOP_REAL, parameter::smooth_fwhm = 30.*coop_SI_arcmin
  COOP_UNKNOWN_STRING,parameter:: color_table = "Rainbow"
  COOP_UNKNOWN_STRING, parameter :: spot_type = "QU"
  COOP_REAL,parameter::r=10.*coop_SI_degree, dr = max(smooth_fwhm/3., r/50.)
  COOP_INT, parameter::n = ceiling(r/dr)

  COOP_UNKNOWN_STRING, parameter :: map_file = "pl353/pl353_iqu.fits"
  COOP_UNKNOWN_STRING, parameter :: spots_file = "spots/pl353_iqu_PmaxSortT_threshold0_fwhm30.txt" 
  COOP_UNKNOWN_STRING, parameter :: imask_file = "predx11/predx11_imask.fits" 
  COOP_UNKNOWN_STRING, parameter :: polmask_file =  "ffp7/ffp7_union_polmask_2048.fits"
  COOP_UNKNOWN_STRING, parameter :: prefix = "multstacked/mult"
  COOP_STRING fout, caption, fname, caption_raw
  COOP_INT,parameter::mmax = 4
  COOP_INT, parameter::npatches = 15
  real, parameter::uppercut  = 3000.
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

  if(index(spots_file, "threshold0").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu$=0"
  elseif(index(spots_file, "threshold1").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu$=1"
  elseif(index(spots_file, "threshold2").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu$=2"
  elseif(index(spots_file, "threshold3").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu$=3"
  elseif(index(spots_file, "threshold4").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu$=4"
  elseif(index(spots_file, "threshold5").gt.0)then
     caption_raw = trim(caption_raw)//", threshold $\nu$=5"
  endif

  call patch(1)%init(spot_type, n, dr, mmax = mmax)
  do i=2, npatches
     patch(i) = patch(1)
  enddo
  if(do_mask)then
     call map%multstack(npatches, patch, -1.e30, uppercut, spots_file, mask)
  else
     call map%multstack(npatches, patch, -1.e30, uppercut, spots_file)
  endif



  do i=1, npatches
     if(i.gt.1)then
        patch(i)%image = (patch(i-1)%image*patch(i-1)%nstack_raw + patch(i)%image*patch(i)%nstack_raw)/(patch(i-1)%nstack_raw + patch(i)%nstack_raw)
        patch(i)%nstack_raw = patch(i-1)%nstack_raw + patch(i)%nstack_raw
     endif
     caption = trim(coop_num2str(patch(i)%nstack_raw))//" patches on "//trim(caption_raw)
     fname = coop_file_name_of(spots_file)
     select case(spot_type)
     case("QrUr")
        fout = prefix//trim(coop_num2str(i))//"_Qr_on_"//trim(fname)
        call patch(i)%plot(imap = 1, output =trim(fout), label = "$Q_r(\mu K)$", caption = caption, color_table = color_table, headless_vectors = .true.)
        call patch(i)%plot(imap = 2, output = prefix//trim(coop_num2str(i))//"_Ur_on_"//trim(fname), label = "$U_r(\mu K)$", caption = caption, color_table = color_table, headless_vectors = .true.)
     case("QU")
        fout = prefix//trim(coop_num2str(i))//"_Q_on_"//trim(fname)
        call patch(i)%plot(imap = 1, output = trim(fout), label = "$Q(\mu K)$", caption = caption, color_table = color_table, headless_vectors = .true.)
        call patch(i)%plot(imap = 2, output = prefix//trim(coop_num2str(i))//"_U_on_"//trim(fname), label = "$U(\mu K)$", caption = caption, color_table = color_table, headless_vectors = .true.)
     case("T", "E", "B", "I") 
        fout = prefix//trim(coop_num2str(i))//"_"//spot_type//"_on_"//trim(fname)
        call patch(i)%plot(imap = 1, output =trim(fout), label = "$"//spot_type//"(\mu K)$", caption = caption, color_table = color_table, headless_vectors = .true.)
     end select
     write(*,*) "the output file is: "//trim(fout)
  enddo
end program test
