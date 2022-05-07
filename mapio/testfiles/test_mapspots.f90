program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use coop_sphere_mod
#ifdef HAS_HEALPIX
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
#endif

  implicit none
#include "constants.h"
  type(coop_file) fp, std
  COOP_INT pix, pix_nest, list(8), nneigh, nbad
  COOP_REAL theta, phi, rot
  COOP_INT, parameter::nside = 1024
  type(coop_healpix_maps)::map
  call fp%open("comparison/maxima_0sigma_ffp7_nobpm_smica_cmb_15a_1024.dat")
  call std%open("spots/comparison_Tmax_threshold0.txt")
!!$  call map%read("comparison/iqu.fits", nmaps_wanted = 1) 
!!$  call map%convert2nested()
  rot = 0.d0
  nbad = 0
  call coop_random_init()
  do 
     read(fp%unit, *, END=100, ERR=100) pix
!!$     call ring2nest(nside, pix, pix_nest)
!!$     call neighbours_nest(nside, pix_nest, list, nneigh)
     call pix2ang_ring(nside, pix, theta, phi)
     call random_number(rot)
     write(std%unit, "(3G18.9)") theta, phi, 0.d0
  enddo
100  call fp%close
  call std%close()
end program test
