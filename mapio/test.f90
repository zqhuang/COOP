program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_file)::fp
  COOP_INT::l
  COOP_REAL::sigmasq 
  sigmasq = (5.d0*coop_SI_arcmin*coop_sigma_by_fwhm)**2/2.d0
  call fp%open("act16/planck_beam.txt")
  do l = 0, 5000
     write(fp%unit, "(I6, E16.7)") l, exp(-l*(l+1.d0)*sigmasq)
  enddo
  call fp%close()
end program test
