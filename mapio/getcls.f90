program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools

  implicit none
#include "constants.h"

  COOP_STRING::imap_file, polmap_file, output_file
  COOP_INT, parameter::lmax = 3000
  COOP_INT, parameter::smooth_delta_ell = 20
  
  COOP_UNKNOWN_STRING, parameter::imask_file = "planck14/dx11_v2_smica_int_mask_005a_2048.fits"
  COOP_UNKNOWN_STRING, parameter::polmask_file = "planck14/dx11_v2_smica_pol_mask_005a_2048.fits"
  COOP_REAL, parameter::rmin = coop_ln2, rmax = log(dble(lmax))
  COOP_REAL lnl(2:lmax)
  type(coop_healpix_maps)::map, imask, polmask
  type(coop_file)::fp
  integer l, i, j
  COOP_INT, parameter::fit_n = 10
  COOP_REAL::coef(fit_n, 3)

  if(iargc().lt.2)then
     write(*,*) "./GetCl imap polmap"
     write(*,*) " or"
     write(*,*) "./GetCl imap polmap outputfile"
     stop
  endif
  imap_file = trim(coop_InputArgs(1))
  polmap_file = trim(coop_InputArgs(2))
  if(iargc() .ge. 3)then
     output_file = trim(coop_InputArgs(3))
  else
     output_file = coop_str_replace(imap_file, ".fits", "_cls.txt")
  endif
  call map%read(imap_file, nmaps_wanted = 3, spin = (/ 0, 2, 2 /) )
  call map%import(polmap_file, index_start = 2, index_end = 3, spin = (/ 2, 2 /))

  call imask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )
  map%map(:,1) = map%map(:, 1)*imask%map(:,1)
  call polmask%read(polmask_file, nmaps_wanted = 1, spin = (/ 0 /) )  
  map%map(:,2) = map%map(:, 2)*polmask%map(:,1)
  map%map(:,3) = map%map(:, 3)*polmask%map(:,1)

  call map%map2alm(lmax = lmax)
  call map%get_cls()
  map%cl(:,1) = map%cl(:,1) * (imask%npix/sum(dble(imask%map(:,1))))
  map%cl(:,2) = map%cl(:,2) * (polmask%npix/sum(dble(polmask%map(:,1))))
  map%cl(:,3) = map%cl(:,3) * (polmask%npix/sum(dble(polmask%map(:,1))))
  call imask%free
  call polmask%free
  
  do l=2, lmax
     lnl(l) = log(dble(l))
  enddo
  coef = 0.d0
  do i=1, 3
     call coop_chebfit(lmax-1, lnl(2:lmax), log(dble(map%cl(2:lmax, i))), fit_n, rmin, rmax, coef(:, i))
     print*,
     do j=1, fit_n-1
        print*, coef(j,i) , ", &"
     enddo
     print*, coef(fit_n, i)
     print*
  enddo
  do l=0, map%lmax
     map%cl(l, :) = map%cl(l, :)*(l*(l+1)/coop_2pi*1.e12)
  enddo
  do i=1, 1
     call coop_smooth_data(map%lmax+1, map%cl(0:map%lmax, i), smooth_delta_ell)
  enddo
  call fp%open(output_file, "w")
  do l = 0, map%lmax
     write(fp%unit, "(I5, 6E16.7)") l, map%cl(l,:)
  enddo
  call fp%close()
!!$  write(*,*) "now producing the lmax filtered map"
!!$  call map%alm2map()
!!$  write(*,*) "nmaps = map%nmaps", map%nmaps, size(map%map, 2)
!!$  write(*,*) "npix = ", map%nside**2*12, map%npix, size(map%map, 1)
!!$  write(*,*) "spin = ", map%spin
!!$  write(*,*) trim(coop_str_replace(output_file, ".txt", ""))//"_I.fits"
!!$  call map%write(trim(coop_str_replace(output_file, ".txt", ""))//"_I.fits", index_list = (/ 1 /) )
!!$  write(*,*) trim(coop_str_replace(output_file, ".txt", ""))//"_QU.fits"  
!!$  call map%write(trim(coop_str_replace(output_file, ".txt", ""))//"_QU.fits", index_list = (/ 2, 3 /) )
end program test
