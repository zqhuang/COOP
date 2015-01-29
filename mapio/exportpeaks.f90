program Exp_spots
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::use_mask = .false.  
  logical::do_max = .true.
  COOP_STRING::peak_name = "$P_T$"
  COOP_STRING::orient_name = "$(Q_T, U_T)$"
  COOP_STRING::map_file = "simu/simu_fullsky_015a_TQTUT_fwhm15.fits" ! "planck14/dx11_v2_smica_int_cmb_010a_1024.fits"
  COOP_STRING::imask_file = "planck14/dx11_v2_common_int_mask_010a_1024.fits"
  COOP_STRING::polmask_file = "planck14/dx11_v2_common_pol_mask_010a_1024.fits"
  COOP_STRING::mask_file_force_to_use = ""
  type(coop_stacking_options)::sto
  type(coop_healpix_maps)::hgm, mask
  COOP_STRING::output = "peaks/simu_imax_nu0"
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
  if(use_mask)then
     if( sto%mask_int .and. .not. sto%mask_pol)then
        call mask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )
     elseif( sto%mask_pol .and. .not. sto%mask_int)then
        call mask%read(polmask_file, nmaps_wanted = 1, spin = (/ 0 /) )
     elseif(trim(mask_file_force_to_use).ne."")then
        call mask%read(mask_file_force_to_use, nmaps_wanted = 1, spin = (/ 0 /) )
     else
        stop "For unknown class of peaks you need to specify the mask file explicitly"
     endif
     if(mask%nside .ne. hgm%nside) stop "mask and map must have the same nside"  
     call hgm%get_peaks(sto, mask = mask)
  else
     call hgm%get_peaks(sto)
  endif
  print*, "find "//COOP_STR_OF(sto%peak_pix%n)//" peaks"
  print*, "writing to "//trim(adjustl(output))//".dat"
  call sto%export(trim(adjustl(output))//".dat")
#else
  print*, "You need to install healpix"
#endif  
end program Exp_spots
