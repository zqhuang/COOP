program test
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
  type(coop_file)::fp
  integer i, nside, id, pix
  type(coop_healpix_maps)::map
  COOP_REAL theta, phi, chisq, prob
  call map%init(nside = 4, nmaps=1, spin = (/ 0 /))
  map%map = 0.
  do i=0, 95
     call fp%open("predx11/predx111024_T_4id"//trim(coop_num2str(i))//".log", "r")
     read(fp%unit, *) nside, id, theta, phi, chisq, prob
     map%map(i, 1) = log10(prob+1.d-4)
     theta = coop_pi - theta
     phi = coop_pi + phi
     call ang2pix_ring(map%nside, theta, phi, pix)
     map%map(pix, 1) = map%map(i, 1)
  enddo
  call map%write("assym.fits")
end program test
