module tmp
  use coop_hnn_mod  
contains
  
  function ur_weight(r)
    real*8 r, ur_weight
    real*8,parameter::r_max = 2. * coop_SI_degree
    ur_weight = sin(r*coop_pi/r_max)**2
  end function ur_weight
end module

program main
  use coop_hnn_mod
  use tmp
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask, mapout
  COOP_INT,parameter::lmax=3000
  COOP_REAL::sqrtCls(0:lmax)
  call map%read("sims/noBB.fits", nmaps_wanted=3)
  mapout = map
  call mask%read("sims/mask256.fits", nmaps_wanted=1)  
  call map%qu2qbarubar(map = mapout, iq = 2, iu = 3,  mask=mask, r_deg = 2.d0, weight=ur_weight)
  call mapout%write("sims/QbarUbar_noBB.fits")
end program main
