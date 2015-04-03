program Exp_spots
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::do_max = .true.
  COOP_STRING::peak_name = "$P$"
  COOP_STRING::orient_name = "$(Q, U)$"
  COOP_STRING::map_file = "planck14/smica_iqu.fits" !! "tuhin/dust_TQUL_015a_b30-500_n512.fits"
  COOP_STRING :: imask_file ="planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_STRING:: polmask_file = "planck14/dx11_v2_common_pol_mask_010a_1024.fits" !"planck14/lat30_mask_n1024.fits" !
  type(coop_stacking_options)::sto
  type(coop_file)::fp
  type(coop_healpix_maps)::hgm, mask, pmap
  COOP_INT::npeaks
  COOP_REAL::theta, phi, angle, rmsI
  COOP_INT i, ip
  call hgm%read(map_file)
  call mask%read(polmask_file)  

  rmsI = sqrt(sum(dble(hgm%map(:,1)**2*mask%map(:,1)))/sum(dble(mask%map(:,1))))
  

  call sto%init(do_max, peak_name, orient_name, nmaps = 3)
  sto%threshold_option = 2
!  sto%I_upper_nu = 55.d0/rmsI
  sto%P_lower_nu = 1.d0
  
  call hgm%get_peaks(sto, mask = mask)
  call sto%export("peaks/smica_pnu1.dat")
  print*, "find ", sto%peak_pix%n, " peaks"

!!$  call pmap%init(nside = 256, nmaps = 1, spin = (/ 0 /) )
!!$  call pmap%mark_peaks(sto, 1)
!!$  
!!$  call pmap%write("dust_peaks.fits", index_list = (/ 1 /) )
#else
  print*, "You need to install healpix"
#endif  
end program Exp_spots
