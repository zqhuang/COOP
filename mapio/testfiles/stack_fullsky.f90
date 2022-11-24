program test
  use coop_wrapper_utils
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, rotmap, stacked_map, stacked_mask, rotmask
  COOP_INT::ipix
  COOP_REAL::theta, phi, psi
  call coop_MPI_init()
  call coop_random_init()
  call map%read("simu/simu_i_16_440a_"//trim(coop_inputArgs(1))//".fits") !"lowl/commander_dx11d2_extdata_temp_cmb_n0016_440arc_v1_cr.fits")
  call mask%read("lowl/commander_dx11d2_mask_temp_n0016_likelihood_v1.fits")
  map%map = map%map*mask%map
  stacked_map = map
  stacked_map%map = 0.
  stacked_mask = mask
  stacked_mask%map = 0.
  call map%map2alm()
  rotmap = map
  rotmask = mask
  do ipix = 0, map%npix-1
     if(map%map(ipix,1).gt. 0.)then
        call map%pix2ang(ipix, theta, phi)
        psi = coop_random_unit()*coop_2pi        
        call map%rotate_coor(output = rotmap, zyz_theta = theta, zyz_phi = phi)
        call rotmap%rotate_coor(zyz_psi =psi)
        call mask%rotate_coor(output = rotmask, zyz_theta = theta, zyz_phi = phi)
        call rotmask%rotate_coor(zyz_psi = psi)
        stacked_map%map = stacked_map%map + rotmap%map
        stacked_mask%map = stacked_mask%map + rotmask%map
     endif
  enddo
  where(stacked_mask%map .gt. 0.5)
     stacked_map%map = stacked_map%map/stacked_mask%map
  elsewhere
     stacked_map%map = 0.
  end where
  call stacked_map%write("stacked_fullsky_"//trim(coop_inputArgs(1))//".fits")
  call coop_MPI_finalize()  
end program test
