program shells
  use coop_hnn_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map
  COOP_INT::l
  COOP_REAL::norm

  call map%read("plmaps/Planck_SMICA_IQU_2048.fits", nmaps_wanted=1)
  call map%map2alm(lmax = 5000)
  call map%get_cls()
  norm = map%cl(2, 1)*6
  do l=2, map%lmax
     map%alm(l, :, 1) = map%alm(l, :, 1) * sqrt(norm/(map%cl(l, 1)*l*(l+1.)))
  enddo
  call map%alm2map()
  call map%write("scale_invariant_CMB.fits")
end program shells
