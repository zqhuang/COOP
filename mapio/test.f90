program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  use coop_fitsio_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::hm, mask
  COOP_SINGLE::rms 
  call hm%open("planck15/COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits", nmaps_wanted = 1)
  call hm%udgrade(nside = 1024)
  call hm%write("planck15/cmb_i.fits")
!!$  call mask%generate_latcut_mask(nside = 1024, latitude_deg = 30.d0, depth_deg = 1.d0)
!!$  call mask%write("planck15/mask_lat30.fits")

end program test
