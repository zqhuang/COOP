program Exp_spots
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::do_max = .true.
  COOP_STRING::peak_name = "$T$"
  COOP_STRING::orient_name = "NULL"
  COOP_STRING::map_file = "planck14/dx11_v2_smica_int_cmb_010a_1024.fits"
  COOP_STRING::imask_file = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_STRING::polmask_file = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"
  type(coop_stacking_options)::sto
  type(coop_healpix_maps)::hgm, mask
  COOP_STRING::output = "peaks/smica_imax_nu0"
  COOP_REAL::threshold_I = 0.d0
  COOP_STRING::line
  if(iargc() .ge. 8)then
     map_file = coop_InputArgs(1)
     imask_file = coop_InputArgs(2)
     polmask_file = coop_InputArgs(3)
     line = coop_InputArgs(4)
     read(line, *) do_max
     peak_name = coop_InputArgs(5)
     orient_name = coop_InputArgs(6)
     output = coop_InputArgs(7)
     line = coop_InputArgs(8)
     read(line, *) threshold_I
  endif
  
  call hgm%read(map_file)
  call sto%init(do_max, peak_name, orient_name, nmaps = hgm%nmaps)
  if(do_max)then
     sto%I_lower_nu = threshold_I
  else
     sto%I_upper_nu = -threshold_I
  endif
  if( sto%index_Q .eq. 0 .and. sto%index_U .eq. 0 )then
     call mask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )
  else
     call mask%read(polmask_file, nmaps_wanted = 1, spin = (/ 0 /) )
  endif
  if(mask%nside .ne. hgm%nside) stop "mask and map must have the same nside"  
  call hgm%get_peaks(sto, mask = mask)
  call sto%export(trim(adjustl(output))//".dat")
#else
  print*, "You need to install healpix"
#endif  
end program Exp_spots
