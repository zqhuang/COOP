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
  integer, parameter::n=10
  type(coop_file)::fp
  integer l
  COOP_REAL::scal(2, 10)
  print*, size(scal)
  stop
  call fp%open(coop_inputArgs(1), "r")
  do l=1, n
     read(fp%unit, *) scal(:, l)
  enddo
  call fp%close()
  print*, sum(scal(1, :))/n, sum(scal(2,:))/n

  
end program test
