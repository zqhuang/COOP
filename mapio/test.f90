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
  call hgm%init(nside = 256, nmaps = 1, spin = (/ 0 /) )
  call hgm%mark_spots("spots/dust_TQUL_015a_b30-500_n512_fwhm15_TQUL_threshold0pt5_2ndthreshold0pt6.txt", imap = 1, width = 30.d0, length = 180.d0)
  call hgm%write("peaks.fits")
  
end program stackth
