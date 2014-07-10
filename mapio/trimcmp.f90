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
  type(coop_healpix_maps)::map, imask, polmask, mapo, poltrim
  integer l, m, il, i, itest
  integer, parameter::lmax = 1000
  type(coop_file)::fp
  COOP_REAL erms
  call map%read("inps/sim2_iqu_nside512.fits", spin = (/ 0, 2 ,2 /) )
  call imask%read("inps/predx11_imask_nside512.fits")
  call polmask%read("inps/predx11_polmask_nside512.fits")
  poltrim = polmask
  call poltrim%trim_mask(real(coop_SI_degree/2.))
  call map%iqu2teb()

  map%map(:, 1) = map%map(:, 1) * imask%map(:, 1)
  map%map(:, 2) = map%map(:, 2) * poltrim%map(:, 1)
  map%map(:, 3) = map%map(:, 3) * poltrim%map(:, 1)
  call map%write("teb_origin.fits")

  call mapo%read("inps/sim2_iqu_nside512.fits", spin = (/ 0, 2 , 2/) )
  call coop_healpix_diffuse_into_mask(mapo, imask, 3.d0*coop_SI_arcmin)
  call coop_healpix_diffuse_into_mask(mapo, polmask, 3.d0*coop_SI_arcmin, .true.)
  call mapo%iqu2teb()
  mapo%map(:, 1) = mapo%map(:, 1) * imask%map(:, 1)
  mapo%map(:, 2) = mapo%map(:, 2) * poltrim%map(:, 1)
  mapo%map(:, 3) = mapo%map(:, 3) * poltrim%map(:, 1)
  call mapo%write("teb_diffuse.fits")
!  call mapo%write("EB_reconstructed.fits", (/2, 3/))
  mapo%map = mapo%map - map%map
  call mapo%write("diff_diffuse.fits")
  print*,"diffuse: ", sqrt(sum(mapo%map(:,2)**2)/mapo%npix), sqrt(sum(mapo%map(:,3)**2)/mapo%npix)


  call mapo%read("inps/sim2_iqu_nside512.fits", spin = (/ 0, 2 , 2/) )
  mapo%map(:, 1) = mapo%map(:, 1) * imask%map(:, 1)
  mapo%map(:, 2) = mapo%map(:, 2) * polmask%map(:, 1)
  mapo%map(:, 3) = mapo%map(:, 3) * polmask%map(:, 1)
  call mapo%iqu2teb()
  mapo%map(:, 1) = mapo%map(:, 1) * imask%map(:, 1)
  mapo%map(:, 2) = mapo%map(:, 2) * poltrim%map(:, 1)
  mapo%map(:, 3) = mapo%map(:, 3) * poltrim%map(:, 1)
  call mapo%write("teb_pseudo.fits")
  mapo%map = mapo%map - map%map
  call mapo%write("diff_pseudo.fits")
  print*,"pesudo: ", sqrt(sum(mapo%map(:,2)**2)/mapo%npix), sqrt(sum(mapo%map(:,3)**2)/mapo%npix)

  call mapo%read("inps/sim2_iqu_nside512_inp0400.fits", spin = (/ 0, 2 , 2/) )
  call mapo%iqu2teb()
  mapo%map(:, 1) = mapo%map(:, 1) * imask%map(:, 1)
  mapo%map(:, 2) = mapo%map(:, 2) * poltrim%map(:, 1)
  mapo%map(:, 3) = mapo%map(:, 3) * poltrim%map(:, 1)
  call mapo%write("teb_inp.fits")
  mapo%map = mapo%map - map%map
  call mapo%write("diff.fits")
  print*,"realization: ", sqrt(sum(mapo%map(:,2)**2)/mapo%npix), sqrt(sum(mapo%map(:,3)**2)/mapo%npix)


  call mapo%read("inps/sim2_iqu_nside512_inp_mean1000.fits", spin = (/ 0, 2 , 2/) )
  call mapo%iqu2teb()
  mapo%map(:, 1) = mapo%map(:, 1) * imask%map(:, 1)
  mapo%map(:, 2) = mapo%map(:, 2) * poltrim%map(:, 1)
  mapo%map(:, 3) = mapo%map(:, 3) * poltrim%map(:, 1)
  call mapo%write("teb_mean.fits")
  mapo%map = mapo%map - map%map
  call mapo%write("meandiff.fits")
  print*,"mean", sqrt(sum(mapo%map(:,2)**2)/mapo%npix), sqrt(sum(mapo%map(:,3)**2)/mapo%npix)

end program test
