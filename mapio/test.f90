program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map, dervs
  COOP_INT::l
  call map%read("planck14/dx11_v2_smica_int_cmb_080a_0128.fits", nmaps_wanted = 6, nmaps_to_read = 1)
  call map%smooth(fwhm = 5.d0*coop_SI_degree, index_list = (/ 1 /) )
  call map%get_QULDD()
  call map%get_dervs(1, dervs, alms_done = .true.)
  call map%write("tquldd1.fits")
  map = dervs
  dervs%map(:, 2) = map%map(:, 2) - map%map(:, 4)
  dervs%map(:, 3) = map%map(:, 3)*2.
  dervs%map(:, 4) = map%map(:, 2) + map%map(:, 4)
  call dervs%write("tquldd2.fits")
end program test
