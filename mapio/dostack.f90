program test
  use coop_healpix_mod
  use coop_wrapper_utils
  implicit none

#include "constants.h"

  COOP_STRING :: spot_type = "B"
  COOP_STRING :: map_file = "simu/simurp2_iqu_10arc_n1024_converted_to_TEB_submap003.fits"
  COOP_STRING:: spots_file ="spots/simurp2_iqu_10arc_n1024_converted_to_TEB_submap003_Bmax_threshold0_fwhm15.txt"
  COOP_STRING :: imask_file = "ffp7/ffp7_imask.fits"
  COOP_STRING:: polmask_file ="commander/commander_polmask.fits"
  COOP_STRING::unit = "muK"


  COOP_UNKNOWN_STRING,parameter:: color_table = "Planck"
  COOP_REAL, parameter::smooth_fwhm = 0.*coop_SI_arcmin
  COOP_REAL,parameter::r=2.*coop_SI_degree, dr = max(smooth_fwhm/3., r/50.)
  COOP_INT, parameter::n = ceiling(r/dr)
  COOP_UNKNOWN_STRING, parameter :: prefix = "stacked/"
  COOP_STRING fout,fout2, caption, fname
  COOP_INT,parameter::mmax = 4
  integer i, m
  type(coop_healpix_maps) map, mask
  type(coop_healpix_patch) patch
  logical::do_mask
  type(coop_asy)::fp
  if(iargc() .ge. 5)then
     map_file = trim(adjustl(coop_InputArgs(1)))
     spots_file = trim(adjustl(coop_InputArgs(2)))
     spot_type = trim(adjustl(coop_InputArgs(3)))     
     imask_file = trim(adjustl(coop_InputArgs(4)))
     polmask_file = trim(adjustl(coop_InputArgs(5)))
     unit = trim(adjustl(coop_InputArgs(6)))
  endif


  select case(trim(spot_type))
  case("T", "I", "zeta")
     if(trim(imask_file) .ne. "")then
        do_mask = .true.
        call mask%read(trim(imask_file), nmaps_wanted = 1 )
     else
        do_mask = .false.
     endif
  case("E","B")
     if(trim(polmask_file) .ne. "")then
        call mask%read(trim(polmask_file), nmaps_wanted = 1)
        do_mask = .true.
     else
        do_mask = .false.
     endif
  case( "QrUr", "QU")
     if(index(map_file, "QTUT").eq.0)then
        if(trim(polmask_file) .ne. "")then
           call mask%read(trim(polmask_file), nmaps_wanted = 1)
           do_mask = .true.
        else
           do_mask = .false.
        endif
     else
        if(trim(imask_file) .ne. "")then
           do_mask = .true.
           call mask%read(trim(imask_file), nmaps_wanted = 1 )
        else
           do_mask = .false.
        endif
     endif
  case default
     stop "Unknown spot_type"
  end select
  call map%read(trim(map_file))
  if(do_mask)then
     do i=1, min(map%nmaps, 3)
        map%map(:, i) = map%map(:, i)*mask%map(:, 1)
     enddo
  endif
  if(trim(unit).eq."K")then
     map%map = map%map*1.e6
  endif

  if(index(spots_file, "_Tmax_QTUTOrient_") .gt. 0 .or. index(spots_file, "TQUmax").gt. 0)then
     caption = "$T$ maxima, oriented"
  elseif(index(spots_file, "_Tmax_").gt.0)then
     caption = "$T$ maxima, random orientation"
  elseif(index(spots_file, "_zetamax_").gt.0)then
     caption = "$\zeta$ maxima, random orientation"
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

  call patch%init(trim(spot_type), n, dr, mmax = mmax)
  if(do_mask)then
     call map%stack(patch, trim(spots_file), mask, do_weight = .true.)
  else
     call map%stack(patch, trim(spots_file), do_weight = .true.)
  endif
  patch%caption = trim(coop_num2str(patch%nstack_raw))//" patches on "//trim(caption)
  patch%color_table = color_table
  fname = coop_file_name_of(trim(spots_file))
  fout2 = ""
  select case(trim(spot_type))
  case("QrUr")
     if(index(map_file, "QTUT") .eq. 0)then
        fout = prefix//"Qr_on_"//trim(fname)
        fout2 = prefix//"Ur_on_"//trim(fname)
     else
        patch%label(1) = "$Q_{T,r}(\mu K)$"
        patch%label(2) = "$U_{T,r}(\mu K)$"
        fout = prefix//"QTr_on_"//trim(fname)
        fout2 = prefix//"UTr_on_"//trim(fname)
     endif
  case("QU")
     if(index(map_file, "QTUT") .eq. 0)then
        fout = prefix//"Q_on_"//trim(fname)
        fout2 = prefix//"U_on_"//trim(fname)
     else
        fout = prefix//"QT_on_"//trim(fname)
        fout2 = prefix//"UT_on_"//trim(fname)
        patch%label(1) = "$Q_T(\mu K)$"
        patch%label(2) = "$U_T(\mu K)$"
     endif
  case("T", "E", "B", "I", "zeta") 
     fout = prefix//trim(spot_type)//"_on_"//trim(fname)
  end select
  call patch%plot(imap = 1, output =trim(fout))
  if(trim(fout2).ne."")call patch%plot(imap = 2, output =trim(fout2))
  write(*,*) "the output file is: "//trim(fout)
  call patch%get_all_radial_profiles()
  do m = 0, 4, 2
     call fp%open(trim(coop_file_add_postfix(fout, "_m"//COOP_STR_OF(m))))
     call fp%init(xlabel="$r$", ylabel="radial profile")
     call coop_asy_curve(fp, patch%r, patch%fr(:, m/2, 1))
     call fp%close()
  enddo
end program test
