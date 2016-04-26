program ClFromMap
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use coop_fitsio_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools

  implicit none
#include "constants.h"

  COOP_STRING::fmap1, fmap2, clfile, fmask
  type(coop_healpix_maps)::map1, map2, mask
  type(coop_file)::fp
  logical::has_mask, has_map2
  COOP_INT::l, lmax, lmax_mask, lmin, nmaps, numcls, i
  COOP_REAL,dimension(:,:),allocatable::Cl_Pseudo, Cl
  COOP_REAL,dimension(:,:,:),allocatable::kernel
  COOP_REAL,dimension(:),allocatable::Cl_mask
  COOP_STRING::genre, unit
  if(iargc().lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./MAP2CLS -map1 ... [-nmaps ...] [-map2 ...] [-mask ...] -lmax ... [-lmax_mask ..] [-lmin ...] -out ..."
     write(*,*) "example:"
     write(*,*) "./MAP2CLS -map1 planck_smica_iqu.fits  -lmax 2000 -lmax_mask 200 -mask planck_smica_mask.fits -fields TEB -out mycls.fits"
     stop
  endif
  call coop_get_command_line_argument(key = "map1", arg = fmap1)
  call coop_get_command_line_argument(key = "nmaps", arg = nmaps, default = 0)
  call coop_get_command_line_argument(key = "lmax", arg = lmax)
  call coop_get_command_line_argument(key = "lmin", arg = lmin, default = 0)
  if(lmin .gt. lmax) stop "lmin > lmax?"
  call coop_get_command_line_argument(key = "map2", arg = fmap2, default = "")
  has_map2 = trim(fmap2).ne.""
  call coop_get_command_line_argument(key = "mask", arg = fmask, default = "")
  has_mask = trim(fmask) .ne. ""
  if(has_mask)then
     call coop_get_command_line_argument(key = "lmax_mask", arg = lmax_mask, default = 200)
  endif
  call coop_get_command_line_argument(key = "out", arg = clfile)
  if(nmaps .gt. 0)then
     call map1%read(fmap1, nmaps_wanted = nmaps)
  else
     call map1%read(fmap1)
     nmaps = map1%nmaps
  endif
  genre = ''
  do i = 1, map1%nmaps
     genre(i:i) = map1%fields(i)(1:1)
  enddo
  numcls = nmaps*(nmaps+1)/2

  if(trim(adjustl(map1%units(1))).ne."")then
     unit = trim(adjustl(map1%units(1)))
  else
     unit = "Unknown"
  endif

  if(has_map2)call map2%read(fmap2, nmaps_wanted = nmaps)

  if(has_mask)then
     call mask%read(fmask, nmaps_wanted = 1)
     allocate(cl_mask(0:lmax_mask))
     call mask%map2alm(lmax = lmax_mask)
     cl_mask = mask%cl(0:lmax_mask, 1)
     allocate(kernel(lmin:lmax, lmin:lmax, 4), Cl_pseudo(lmin:lmax, 6), Cl(lmin:lmax, 6))
     Cl = 0.d0
     Cl_pseudo = 0.d0

     call map1%apply_mask(mask)
     call map1%map2alm(lmax = lmax)
     if(has_map2)then
        call map2%apply_mask(mask)
        call map2%map2alm(lmax = lmax)
        call map1%get_cls(map2)
     else
        call map1%get_cls()
     endif
     select case(trim(genre))
     case("I", "T")
        call coop_pseudoCl_get_kernel(lmax_mask, Cl_mask, lmin, lmax, kernel, want_T = .true., want_EB = .false.)
        cl_pseudo(lmin:lmax, coop_TEB_index_TT) = map1%cl(lmin:lmax, 1)
        call coop_pseudoCl2Cl(lmin = lmin, lmax = lmax, Cl_pseudo = Cl_pseudo, kernel = kernel, Cl = Cl, want_T = .true., want_EB = .false.)
        call coop_fits_file_write_cls(lmin = lmin, lmax =lmax, cls = cl(lmin:lmax,  coop_TEB_index_TT:coop_TEB_index_TT ), filename = clfile, genre = genre, spin = (/ 0 /) , unit=unit)
     case("QU")
        call coop_pseudoCl_get_kernel(lmax_mask, Cl_mask, lmin, lmax, kernel, want_T = .false., want_EB = .true.)
        cl_pseudo(lmin:lmax, coop_TEB_index_EE) = map1%cl( lmin:lmax, COOP_MATSYM_INDEX(2, 1, 1) )
        cl_pseudo(lmin:lmax, coop_TEB_index_BB) = map1%cl( lmin:lmax,  COOP_MATSYM_INDEX(2, 2, 2) )
        cl_pseudo(lmin:lmax, coop_TEB_index_EB) = map1%cl( lmin:lmax,  COOP_MATSYM_INDEX(2, 1, 2) )
        call coop_pseudoCl2Cl(lmin = lmin, lmax = lmax, Cl_pseudo = Cl_pseudo, kernel = kernel,  Cl = Cl, want_T = .false., want_EB = .true.)
        call coop_fits_file_write_cls(lmin = lmin, lmax =lmax, cls = cl(lmin:lmax, (/ coop_TEB_index_EE, coop_TEB_index_BB, coop_TEB_index_EB /) ), filename = clfile, genre = genre, spin = (/ 2, 2 /) , unit=unit)
     case("IQU", "TQU")
        call coop_pseudoCl_get_kernel(lmax_mask, Cl_mask, lmin, lmax, kernel, want_T = .true., want_EB = .true.)
        cl_pseudo(lmin:lmax, 1:6) = map1%cl(lmin:lmax, 1:6)
        call coop_pseudoCl2Cl(lmin = lmin, lmax = lmax, Cl_pseudo = Cl_pseudo, kernel = kernel,  Cl = Cl, want_T = .true., want_EB = .true.)
        call coop_fits_file_write_cls(lmin = lmin, lmax =lmax, cls = cl, filename = clfile, genre = genre, spin = (/ 0, 2, 2 /) , unit=unit)
     case default
        stop "In the mask mode MAP2CLS only supports fields = I, QU, or IQU"
     end select

  else
     call map1%map2alm(lmax = lmax)
     if(has_map2)then
        call map2%map2alm(lmax = lmax)
        call map1%get_cls(map2)
     else
        call map1%get_cls()
     endif
     allocate(cl(lmin:lmax, numcls))
     cl = map1%cl(lmin:lmax, 1:numcls)
     write(*,*) trim(genre)
     write(*,*) trim(unit)
     call coop_fits_file_write_cls(lmin = lmin, lmax =lmax, cls = cl, filename = trim(clfile),  genre = genre, spin = map1%spin(1:nmaps), unit=unit)
  endif

end program ClFromMap
