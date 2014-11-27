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
  integer, parameter::lmax=1500
  type(coop_file)::fp
  integer l, il
  COOP_REAL::scal(2, lmax)
  call fp%open("calcls.txt", "r")
  do l=2, 1500
     read(fp%unit, *) il, scal(:, l)
  enddo
  call fp%close()

  
end program test
