program stackth
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
  type(coop_healpix_maps)::hgm, mask, hgm2
  COOP_INT l
  
  call hgm%read("tuhin/dust_iqu_512.fits", nmaps_wanted = 3, spin=(/ 0, 2, 2 /) )
  call mask%read("planck14/lat25_mask_n512.fits")
  hgm%map(:,2) = hgm%map(:,2)*mask%map(:,1)
  hgm%map(:,3) = hgm%map(:,3)*mask%map(:,1)
  call hgm%map2alm( index_list = (/ 2, 3 /) )
  do l=0, hgm%lmax
     hgm%alm(l, :, 2:3) = hgm%alm(l, :, 2:3) * (coop_highpass_filter(40, 60, l)*coop_lowpass_filter(250, 350, l)*coop_gaussian_filter(15.d0, l))
  enddo
  call hgm%alm2map( index_list = (/2, 3/))
  call hgm%write("tuhin/dust_QU_015a_512.fits", index_list = (/ 2 ,3 /) )
  call hgm%iqu2TEB()
  call hgm%write("tuhin/dust_EB_015a_512.fits", index_list=(/2,3/))
  
end program stackth
