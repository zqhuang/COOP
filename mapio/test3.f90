program shells
  use coop_hnn_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map
  call map%generate_latcut_mask(nside = 512, latitude_deg = 35.d0)
  call map%write("sims/mask512_l30.fits")
end program shells
