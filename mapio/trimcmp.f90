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
  type(coop_healpix_maps)::map, imask, polmask, mapo
  integer l, m, il, i
  integer, parameter::lmax = 1000
  type(coop_file)::fp
  COOP_REAL erms
  call map%read("inps/simu_iqu_nside512_inp_mean0400.fits", spin = (/ 0, 2 ,2 /) )
  call imask%read("inps/predx11_imask_nside512.fits")
  call polmask%read("inps/predx11_polmask_nside512.fits")
!  call polmask%trim_mask(coop_SI_degree)
  call map%iqu2teb()
  
  map%map(:, 1) = map%map(:, 1) * imask%map(:, 1)
  map%map(:, 2) = map%map(:, 2) * polmask%map(:, 1)
  map%map(:, 3) = map%map(:, 3) * polmask%map(:, 1)
  call map%write("teb_inp.fits")

  call mapo%read("inps/simu_iqu_nside512.fits", spin = (/ 0, 2 , 2/) )
  call mapo%map2alm()
  call mapo%iqu2teb()
  mapo%map(:, 1) = mapo%map(:, 1) * imask%map(:, 1)
  mapo%map(:, 2) = mapo%map(:, 2) * polmask%map(:, 1)
  mapo%map(:, 3) = mapo%map(:, 3) * polmask%map(:, 1)
  call mapo%write("teb_origin.fits")
  erms = sqrt(sum(mapo%map(:,2)**2)/sum(polmask%map))
  map%map = map%map - mapo%map
  call map%write("diff2.fits")
  map%map = map%map/erms
  call map%write("reldiff2.fits")

end program test
