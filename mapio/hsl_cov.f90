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

  COOP_UNKNOWN_STRING, parameter::prefix="T_on_Tmax_"
  COOP_INT, parameter::nside = 4
  COOP_INT, parameter::nfiles = nside**2*6
  COOP_INT, parameter::nsims = 1000
  COOP_INT, parameter::npix = 30
  COOP_INT i, j, k1, k2
  COOP_REAL::diff(npix, nsims), cov(npix, npix), ddf(npix)
  type(coop_file)::fp
  do i=0, nfiles-1
     call fp%open(prefix//"fr_"//COOP_STR_OF(nside)//"_"//COOP_STR_OF(i)//".txt", "r")
     do j=1, nsims
        read(fp%unit, *) diff(:, j)
     enddo
     call fp%close()
     
  enddo





end program test
