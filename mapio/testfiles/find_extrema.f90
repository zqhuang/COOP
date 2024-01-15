program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  use udgrade_nr
  use coord_v_convert,only:coordsys2euler_zyz
#endif  
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::hmap, lmap, mmap, mask, mout
  COOP_INT::hottest_pix, hotpix(20), listpix(0:12*256**2-1), nlist(8), nn, locs(3), locs2(4), nhot
  COOP_REAL::tn(8), cut
  COOP_INT::i, j, ipix
  call lmap%read("lowl/commander_0002.fits")
  call hmap%read("lowl/commander_dx11d2_extdata_temp_cmb_n0256_60arc_v1_cr_smoothed_fwhm600arcmin.fits")
  call mask%read("lowl/commander_dx11d2_mask_temp_n0256_likelihood_v1.fits")
  call lmap%convert2nested()
  call hmap%convert2nested()  
  call mask%convert2nested()
  call mout%init(nside = hmap%nside, nmaps = 1, genre = "MASK", nested=.true.)
  mout%map = 0.
  hotpix(1) = coop_maxloc(lmap%map(:,1))-1
  call neighbours_nest(lmap%nside, hotpix(1), nlist, nn)
  tn(1:nn) = lmap%map(nlist(1:nn), 1)
  call coop_find_maxlocs(tn(1:nn), locs)
  hotpix(2:4) = nlist(locs)
  call neighbours_nest(lmap%nside, hotpix(2), nlist, nn)  
  tn(1:nn) = lmap%map(nlist(1:nn), 1)
  call coop_find_maxlocs(tn(1:nn), locs2)
  hotpix(5:6) = nlist(locs2(1:2))
  lmap%map = 0.
  lmap%map(hotpix(1:6),1) = (/ 1., 2., 3., 4., 5., 6. /)
  call lmap%write("hot.fits")
  cut =0.! minval(lmap%map(hotpix,1))
  do i=1, 6
     do j = hotpix(i)*(hmap%npix/lmap%npix), (hotpix(i)+1)*(hmap%npix/lmap%npix)-1
        if(hmap%map(j, 1) .gt. cut)then
           mout%map(j,1) = 1.
        endif
     enddo
  enddo
  mout%map = mout%map*mask%map
  call mout%write("lowl/mask_hot_bar_n0256.fits")
end program test
