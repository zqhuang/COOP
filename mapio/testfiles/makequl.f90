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
  COOP_INT,parameter::llow = 50, lhigh =400
  call hgm%read("tuhin/dust_iqu_512.fits", nmaps_wanted = 4, nmaps_to_read = 1, spin=(/ 0, 2, 2, 0 /) )
  call mask%read("planck14/lat30_mask_s512.fits")
  hgm%map(:,1) = hgm%map(:,1)*mask%map(:,1)/287.d0 !!do mask, convert to K
  call hgm%map2alm( index_list = (/ 1 /) )
  hgm%alm(:, :, 3) = 0.  
  do l=0, hgm%lmax
     hgm%alm(l, :, 1) = hgm%alm(l, :, 1) * (coop_highpass_filter(llow-10, llow+10, l)*coop_lowpass_filter(lhigh-50, lhigh+50, l)*coop_gaussian_filter(15.d0, l))
     hgm%alm(l, :, 2) = hgm%alm(l,:,1)*l*(l+1.d0)
     hgm%alm(l, :, 4) = hgm%alm(l,:,2)
  enddo
  call hgm%alm2map()
  call hgm%write("tuhin/dust_TQUL_015a_b"//COOP_STR_OF(llow)//"-"//COOP_STR_OF(lhigh)//"_n512.fits")
  
end program stackth
