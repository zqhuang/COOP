program Exp_spots
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::do_max = .true.
  COOP_STRING::peak_name = "RANDOM"
  COOP_STRING::orient_name = "$(\nabla^2 Q_T, \nabla^2 U_T)$"
  COOP_STRING::map_file =  "tuhin/dust_TQUL_015a_b30-500_n512.fits"
  COOP_STRING :: imask_file = "planck14/lat30_mask_n512.fits" ! "planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_STRING:: polmask_file = "planck14/lat30_mask_n512.fits" ! "planck14/dx11_v2_common_pol_mask_010a_1024.fits" !
  type(coop_stacking_options)::sto
  type(coop_file)::fp
  type(coop_healpix_maps)::hgm, mask, pmap
  COOP_INT::npeaks
  COOP_REAL::theta, phi, angle, rmsI
  COOP_INT i, ip
  call hgm%read(map_file)
  call mask%read(polmask_file)  

  rmsI = sqrt(sum(dble(hgm%map(:,1)**2*mask%map(:,1)))/sum(dble(mask%map(:,1))))

  call sto%init(do_max, peak_name, orient_name, nmaps = 4)
  sto%threshold_option = 7
  sto%I_lower_nu = 6.d-5/rmsI
  sto%P2byL2_lower = 0.**2
  sto%P2byL2_upper = 0.2d0**2


  call hgm%get_peaks(sto, mask = mask)
  print*, "find "//COOP_STR_OF(sto%peak_pix%n)//" peaks"
  call sto%export("peaks/dust_fp_e0to0pt2.dat")  
!!$  call pmap%init(nside = 256, nmaps = 1, spin = (/ 0 /) )
!!$  call pmap%mark_peaks(sto, 1)
!!$  
!!$  call pmap%write("dust_peaks.fits", index_list = (/ 1 /) )
#else
  print*, "You need to install healpix"
#endif  
end program Exp_spots
