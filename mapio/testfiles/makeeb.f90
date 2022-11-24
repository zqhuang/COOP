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
  COOP_INT,parameter::llow = 50
  COOP_INT,parameter::lhigh = 400
  call hgm%read("tuhin/dust_iqu_512.fits", nmaps_wanted = 4, nmaps_to_read = 3, spin=(/ 0, 2, 2, 0 /) )
  call mask%read("planck14/lat30_mask_s512.fits")
  hgm%map(:,2) = hgm%map(:,2)*mask%map(:,1)
  hgm%map(:,3) = hgm%map(:,3)*mask%map(:,1)  
  call hgm%map2alm( index_list = (/ 2, 3 /) )
  do l=0, hgm%lmax
     hgm%alm(l, :, 2:3) = hgm%alm(l, :, 2:3) * (coop_highpass_filter(llow-10, llow+10, l)*coop_lowpass_filter(lhigh-50, lhigh+50, l)*coop_gaussian_filter(15.d0, l))
  enddo
  call hgm%alm2map(index_list = (/ 2, 3 /) )
  call hgm%write("tuhin/dust_QU_015a_b"//COOP_STR_OF(llow)//"-"//COOP_STR_OF(lhigh)//"_n512.fits", index_list = (/ 2, 3/) )
  hgm%spin = 0
  call hgm%alm2map(index_list = (/ 2, 3 /) )
  call hgm%write("tuhin/dust_E_015a_b"//COOP_STR_OF(llow)//"-"//COOP_STR_OF(lhigh)//"_n512.fits", index_list = (/ 2 /) )
  call hgm%write("tuhin/dust_B_015a_b"//COOP_STR_OF(llow)//"-"//COOP_STR_OF(lhigh)//"_n512.fits", index_list = (/ 3 /) )  

  
end program stackth
