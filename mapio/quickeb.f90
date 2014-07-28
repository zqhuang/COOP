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
  type(coop_healpix_maps)::map, imask, polmask, mapo, poltrim, tmp
  integer l, itest
  COOP_REAL rms
  integer, parameter::lmax = 1000
  type(coop_file)::fp
  call map%read("ffp7/ffp7_cmb_plus_noise_0001_15a_1024.fits", spin = (/ 0, 2 ,2 /) )
  tmp = map
  call imask%read("ffp7/ffp7_imask1024.fits")
  call polmask%read("ffp7/ffp7_union_polmask_1024.fits")
  poltrim = polmask
  call poltrim%trim_mask(real(coop_SI_degree))
  call map%iqu2teb()

  call map%mask(imask, polmask)
!  call map%map2alm()

  rms =sqrt(sum( map%map(:,2)**2)/map%npix)

  do itest = 0, 12
     mapo = tmp
     call mapo%mask(imask, polmask)
     if(itest.gt.0)call coop_healpix_diffuse_into_mask(mapo, polmask, coop_SI_arcmin*itest*20, .true.)
     call mapo%iqu2teb()
     call mapo%mask(imask, polmask)
     print*, itest*20,sqrt(sum((mapo%map(:,2)-map%map(:,2))**2)/mapo%npix)/rms
  enddo

!!$  call mapo%map2alm()
!!$
!!$  call fp%open("pseudoCls_laplacian.txt")
!!$  do l = 2, 2000
!!$     write(fp%unit, "(I5, 2E16.7)") l, map%cl(l, coop_healpix_index_EE), mapo%cl(l, coop_healpix_index_EE)
!!$  enddo
!!$  call fp%close()

end program test
