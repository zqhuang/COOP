program Exp_spots
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::do_max = .true.
  COOP_STRING::peak_name = "$T$"
  COOP_STRING::orient_name = "$\nabla^2$"
  COOP_STRING::map_file = "tuhin/dust_TQUL_015a_b30-500_n512.fits"
  COOP_STRING :: imask_file = "planck14/lat30_mask_n512.fits"
  COOP_STRING:: polmask_file = "planck14/lat30_mask_n512.fits"
  type(coop_stacking_options)::sto
  type(coop_file)::fp
  type(coop_healpix_maps)::hgm, mask, pmap
  COOP_INT::npeaks
  COOP_REAL::theta, phi, angle
  COOP_INT i, ip
  call sto%init(do_max, peak_name, orient_name, nmaps = 4)
  sto%threshold_option = 7
  sto%I_lower_nu = 0.5d0
  sto%P2byL2_lower = 0.6d0**2
  call mask%read(imask_file)
  call hgm%read(map_file)
  call hgm%get_peaks(sto, mask = mask)
  call sto%export_ang("catlog.txt")
  print*, "find ", sto%peak_pix%n, " peaks"
  call pmap%init(nside = 256, nmaps = 1, spin = (/ 0 /) )
  call pmap%mark_peaks(sto, 1)
  
  call pmap%write("dust_peaks.fits", index_list = (/ 1 /) )
#else
  print*, "You need to install healpix"
#endif  
end program Exp_spots
