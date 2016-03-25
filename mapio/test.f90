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
  call hm%open("planck15/cmb_i.fits")
  call mask%generate_latcut_mask(nside = 1024, latitude_deg = 30.d0, depth_deg = 1.d0)


end program test
