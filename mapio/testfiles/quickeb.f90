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
  type(coop_healpix_maps)::map, imask, polmask, mapo, poltrim, tmp, fullmap
  integer l, itest
  COOP_REAL rms
  integer, parameter::lmax = 1000
  type(coop_file)::fp
  call map%read("ffp7/ffp7_cmb_plus_noise_0001_15a_1024.fits", spin = (/ 0, 2 ,2 /) )
  map%map = map%map
  mapo = map
  fullmap = map

  call imask%read("ffp7/ffp7_imask1024.fits")
  call polmask%read("ffp7/ffp7_union_polmask_1024.fits")
  poltrim = polmask
!  call poltrim%trim_mask(real(coop_SI_degree))
  call map%iqu2teb()

  call map%mask(imask, poltrim)
  call map%map2alm()
  call fullmap%map2alm

  rms =sqrt(sum( map%map(:,2)**2)/map%npix)

  call mapo%mask(imask, polmask)
  call coop_healpix_diffuse_into_mask(mapo, polmask, coop_SI_arcmin*30, .true.)
  call mapo%iqu2teb()
  call mapo%write("ffp7_nobpm_smica_cmbPlusnoise_0001_15a_1024_EB.fits", (/2, 3/) )
  call mapo%mask(imask, poltrim)
  call mapo%map2alm()
  tmp = mapo
  tmp%map = tmp%map - map%map
  call tmp%map2alm()
  call fp%open("pseudoCls.txt")
  do l = 2, 1500
     write(fp%unit, "(I5, 10E16.7)") l,  tmp%cl(l, coop_healpix_index_EE)/map%cl(l, coop_healpix_index_EE)
  enddo
  call fp%close()

end program test
