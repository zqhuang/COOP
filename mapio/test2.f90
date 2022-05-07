program shells
  use coop_hnn_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, mask
  integer i
  call map%read("plmaps/DustFilaments_TQU_NS2048_Nfil180p5M_CMBS4BICEP_353p0GHz.fits", nmaps_wanted=3)

  call mask%read("sims/mask_n2048_l30.fits")
  do i=1, map%nmaps
     map%map(:, i) = map%map(:, i) * mask%map(:, 1)
  enddo
  call map%smooth(fwhm = 15.d0*coop_SI_arcmin, l_lower = 80, delta_l = 10)
  call map%write("sims/DustFil_353_fwhm15_lowl80.fits")
end program shells
