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
  COOP_INT l
  do l=2, 2500
     print*, l, coop_Planck_TNoise(l)*(l*(l+1)/coop_2pi*2.726e6**2), coop_Planck_ENoise(l)*(l*(l+1)/coop_2pi*2.726e6**2), coop_Planck_BNoise(l)*(l*(l+1)/coop_2pi*2.726e6**2)
  enddo
end program stackth
