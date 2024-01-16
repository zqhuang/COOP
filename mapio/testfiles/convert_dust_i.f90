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
  call hgm%read("tuhin/dust_i_512.fits", nmaps_wanted = 3, nmaps_to_read = 1, spin=(/ 0, 2, 2 /) )
!  call mask%read("planck14/dx11_v2_smica_int_mask_010a_1024.fits")
  !  hgm%map(:,1) = hgm%map(:,1)*mask%map(:,1)
  hgm%map(:,1) = hgm%map(:,1)/287.d0 !!convert to K
  call hgm%map2alm( index_list = (/ 1 /) )
  do l=0, hgm%lmax
     hgm%alm(l, :, 1) = hgm%alm(l, :, 1) * (coop_highpass_filter(40, 60, l)*coop_lowpass_filter(250, 350, l)*coop_gaussian_filter(15.d0, l))
  enddo
  call hgm%alm2map(index_list = (/ 1 /) )
  call hgm%iqu2TQTUT(idone = .true.)  
  call hgm%write("tuhin/dust_TQTUT_015a_512.fits")
  hgm2 = hgm
  
  do l=0, min(hgm%lmax, 550)
     hgm%alm(l, :, 1) = hgm%alm(l, :, 1) * l * (l+1.d0)
  enddo
  call hgm%alm2map(index_list = (/ 1 /) )  
  call hgm%iqu2TQTUT(idone = .true.)
  call hgm%write("tuhin/dust_lap_TQTUT_015a_512.fits")
  hgm%map(:,1) = hgm2%map(:,1)
  call hgm%write("tuhin/dust_TQTUT_lap_015a_512.fits")
  
end program stackth
