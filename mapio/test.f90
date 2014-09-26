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
  integer l
  type(coop_healpix_maps)::map
  call map%read("planck/smica_inp_cmb.fits")
  call map%map2alm()
  do l = 2, map%lmax
     map%alm(l,:,1) = map%alm(l,:,1)*exp(-((l-25.)/5.)**2)
  enddo
  call map%alm2map()
  call map%write("mapdeficit.fits")
  
end program test
