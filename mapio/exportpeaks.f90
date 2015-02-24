program Exp_spots
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::use_mask = .true.
  logical::do_max = .false.
  logical::remove_mono = .false.  
  COOP_STRING::peak_name = "$T$"
  COOP_STRING::orient_name = "NULL"
  COOP_STRING::map_file =  "act15/act15_i_hp_230_270_smoothed_fwhm5arcmin.fits"
  COOP_STRING::imask_file = "act15/act15_imask.fits"
  COOP_STRING::polmask_file = "act15/act15_polmask.fits"
  COOP_STRING::mask_file_force_to_use = ""
  
  type(coop_stacking_options)::sto
  type(coop_healpix_maps)::hgm, mask
  COOP_STRING::output = "peaks/act_imin_5a"
  COOP_REAL::threshold = 0.
  COOP_STRING::line
  COOP_INT::i
  
  if(iargc() .ge. 8)then
     use_mask = .true.
     map_file = coop_InputArgs(1)
     imask_file = coop_InputArgs(2)
     polmask_file = coop_InputArgs(3)
     line = coop_InputArgs(4)
     read(line, *) do_max
     peak_name = coop_InputArgs(5)
     orient_name = coop_InputArgs(6)
     output = coop_InputArgs(7)
     line = coop_InputArgs(8)
     read(line, *) threshold
  else
     coop_healpix_mask_tol = 0.
     if(.not. use_mask)then
        write(*,*) "Warning: not using the mask"
     endif
  endif
  
  call hgm%read(map_file)

  call sto%init(do_max, peak_name, orient_name, nmaps = hgm%nmaps)
  select case(trim(coop_str_numalpha(peak_name)))
  case("T", "I", "E", "B", "zeta", "Z")
     if(do_max)then
        sto%I_lower_nu = threshold
     else
        sto%I_upper_nu = -threshold
     endif
  case("P", "PT", "Pzeta", "PZ")
     if(do_max)then
        sto%P_lower_nu = max(threshold, 0.d0)
     else
        if(threshold .le. 0.d0) stop "For Pmin you must use positive threshold"
        sto%P_upper_nu = threshold
     endif
  case default
     if(abs(threshold).lt. coop_stacking_max_threshold)then
        write(*,*) "cannot apply threshold for unknown peak name: "//trim(peak_name)
        stop
     endif
  end select
  if(use_mask)then
     if( sto%mask_int .and. .not. sto%mask_pol)then
        call mask%read(imask_file, nmaps_wanted = 1, spin = (/ 0 /) )
     elseif( sto%mask_pol .and. .not. sto%mask_int)then
        call mask%read(polmask_file, nmaps_wanted = 1, spin = (/ 0 /) )
     elseif(trim(mask_file_force_to_use).ne."")then
        call mask%read(mask_file_force_to_use, nmaps_wanted = 1, spin = (/ 0 /) )
        write(*,*) "Using forced mask: "//trim(mask_file_force_to_use)
     else
        stop "For unknown class of peaks you need to specify the mask file explicitly"
     endif
     if(mask%nside .ne. hgm%nside) stop "mask and map must have the same nside"
     if(remove_mono)then
        do i = 1, hgm%nmaps
           hgm%map(:,i) =  hgm%map(:,i)* mask%map(:,1)
           hgm%map(:,i) =  hgm%map(:,i) - sum(dble( hgm%map(:,i)))/sum(dble(mask%map(:,1)))
        enddo
     endif
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
