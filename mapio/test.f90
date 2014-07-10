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
  type(coop_healpix_maps)::map, imask, polmask !imap, polmap, inoise, polnoise, iqumap
  call map%read("pl353/353GHz_ful.fits")
  call imask%read("predx11/predx11_imask.fits")
!  call polmask%read("ffp7/ffp7_union_polmask_2048.fits")
  call polmask%read("predx11/predx11_polmask.fits")
  map%map(:,1) = map%map(:,1)*imask%map(:,1)
  map%map(:,2) = map%map(:,2)*polmask%map(:,1)
  map%map(:,3) = map%map(:,3)*polmask%map(:,1)
  map%map = map%map*1.d6
  call map%write("pl353/pl353_iqu.fits")
!!$  call iqumap%init(nside=1024, nmaps = 3, spin = (/ 0, 2, 2 /) )
!!$  call imap%read("pl353/353GHz_ful.fits")
!!$  call inoise%read("ffp7/ffp7_smica_noise_1024_I_smoothed_fwhm15arcmin.fits")
!!$  iqumap%map(:,1:1) = imap%map + inoise%map
!!$  call polmap%read("ffp7/ffp7_nobpm_smica_cmb_mc_0001_15a_1024_QU.fits")
!!$  call polnoise%read("ffp7/ffp7_nobpm_smica_noise_mc_0001_15a_1024_QU.fits")
!!$  iqumap%map(:,2:3) = polmap%map + polnoise%map
!!$  call iqumap%write("ffp7/ffp7_cmb_plus_noise_0001_15a_1024.fits")
end program test
