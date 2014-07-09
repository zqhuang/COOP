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
  call imask%read("inps/predx11_imask_nside512.fits")
  call polmask%read("inps/predx11_polmask_nside512.fits")
!!  call polmask%trim_mask(real(5.*coop_SI_degree))

!!  call map%read("inps/simu_iqu_nside512_inp_mean0400.fits", spin = (/ 0, 2, 2 /) )
  call mapo%read("inps/simu_iqu_nside512.fits", spin = (/ 0, 2 , 2/) )
  call mapo%iqu2teb()
  call mapo%write("test.fits")
  call map%read("test.fits", spin=(/ 0, 0, 0 /) )
  call map%teb2iqu()
  call mapo%read("inps/simu_iqu_nside512.fits", spin = (/ 0, 2 , 2/) )

  map%map = map%map - mapo%map


  map%map(:, 1) = map%map(:, 1) * imask%map(:, 1)
  map%map(:, 2) = map%map(:, 2) * polmask%map(:, 1)
  map%map(:, 3) = map%map(:, 3) * polmask%map(:, 1)
  print*, maxval(abs(map%map(:,1)))
  print*, maxval(abs(map%map(:,2)))
  print*, maxval(abs(map%map(:,3)))
  call map%write("diffqu.fits")
  
!!$  call map%map2alm()
!!$
!!$  call mapo%read("inps/sim2_iqu_nside512.fits", spin = (/ 0, 2 , 2/) )
!!$  call mapo%iqu2lapteb()
!!$  mapo%map(:, 1) = mapo%map(:, 1) * imask%map(:, 1)
!!$  mapo%map(:, 2) = mapo%map(:, 2) * polmask%map(:, 1)
!!$  mapo%map(:, 3) = mapo%map(:, 3) * polmask%map(:, 1)
!!$  call mapo%map2alm()
!!$
!!$  call fp%open("pseudo_cls.txt", "w")
!!$  do l=2, 1000
!!$     write(fp%unit, "(I5, 2E16.7)") l, map%cl(l, coop_healpix_index_EE),  mapo%cl(l, coop_healpix_index_EE)
!!$  enddo
!!$  call fp%close()

end program test
