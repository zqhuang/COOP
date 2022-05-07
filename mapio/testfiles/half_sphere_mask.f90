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

  COOP_UNKNOWN_STRING, parameter:: input_mask = "predx11/predx11_imask.fits" ! "ffp7/ffp7_smica_harmonic_mask_I.fits"
  COOP_REAL, parameter:: axis_l =226.d0
  COOP_REAL, parameter:: axis_b = -17.d0
  type(coop_healpix_maps)::hm, hmnorth
  integer pix, nlist
  integer, dimension(:),allocatable::listpix
  COOP_REAL vec(3), theta, phi
  call hm%read(input_mask)
  allocate(listpix(0:hm%npix-1))
  call coop_healpix_lb2ang(axis_l, axis_b, theta, phi)
  call ang2vec(theta, phi, vec)
  call query_disc(hm%nside, vec, coop_pio2, listpix, nlist, nest = 0, inclusive = 0)
  hmnorth = hm
  hm%map = 0.
  hm%map(listpix(0:nlist-1), :) = hmnorth%map(listpix(0:nlist-1), :)
  hmnorth%map(listpix(0:nlist-1), :) = 0.
  call hm%write(trim(coop_file_add_postfix(input_mask, "_south")))
  call hmnorth%write(trim(coop_file_add_postfix(input_mask, "_north")))
  print*, sum(hm%map), sum(hmnorth%map)
end program test
