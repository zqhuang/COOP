program Exp_spots
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::remove_mono = .true.  
  logical::hot, cold
  COOP_STRING::map_file, imask_file, polmask_file, mask_file, peak_name, orient_name, output
  type(coop_stacking_options)::sto
  type(coop_healpix_maps)::hgm, mask
  COOP_REAL::threshold
  COOP_INT::i, nside_scan


  if(iargc().le.0)then
     write(*,"(A)") "----------------------------------------------------------"     
     write(*,"(A)") "Syntax:"
     write(*,"(A)") "./GetPeaks -map MAP_FILE -out OUTPUT_FILE"
     write(*,"(A)") "----------------------------------------------------------"
     write(*,"(A)") "other options are:"     
     write(*,"(A)") "-hot [T|F]"
     write(*,"(A)") "-cold [F|T]"
     write(*,"(A)") "-peak [SADDLE|COL|RANDOM|T|E|B|P_T|P|\zeta|P_\zeta]"
     write(*,"(A)") "-orient [RANDOM|(Q_T, U_T)|(Q, U)|(Q_{\nabla^2T},U_{\nabla^2T})]"
     write(*,"(A)") "-mask [MASK_FILE|NONE|'']"
     write(*,"(A)") "-nu THRESHOLD_VALUE"
     write(*,"(A)") "-nss NSIDE_SCAN"          
     write(*,"(A)") "-imask [INTENSITY_MASK_FILE] (automatically chosen when mask = '')"
     write(*,"(A)") "-polmask [POL_MASK_FILE] (automatically chosen when mask = '')"
     write(*,"(A)") "----------------------------------------------------------"     
     stop
  endif
  call coop_get_command_line_argument(key = 'map', arg = map_file)
  call coop_get_command_line_argument(key = 'hot', arg = hot, default =.true.)
  call coop_get_command_line_argument(key = 'cold', arg = cold, default =.not. hot)
  if(.not. hot .and. .not. cold) stop "Error: you cannot set both hot and cold to be false"

  call coop_get_command_line_argument(key = 'mask', arg = mask_file, default = 'NONE')
  call coop_get_command_line_argument(key = 'imask', arg = imask_file, default = '')
  call coop_get_command_line_argument(key = 'polmask', arg = polmask_file, default = '')      
  call coop_get_command_line_argument(key = 'peak', arg = peak_name, default = 'T')
  call coop_get_command_line_argument(key = 'orient', arg = orient_name, default = 'RANDOM')
  call coop_get_command_line_argument(key = 'nu', arg = threshold, default = 0.d0)
  if((index(orient_name, coop_backslash).ne.0 .or. index(orient_name, "_") .ne. 0 .or. index(orient_name, "^").ne.0) .and. index(orient_name, "$").eq.0)then
     orient_name = "$"//trim(adjustl(orient_name))//"$"
  endif
  call coop_get_command_line_argument(key = 'out', arg = output)
  call coop_get_command_line_argument(key = 'nss', arg = nside_scan, default = 0)
  coop_healpix_mask_tol = 0.
  call hgm%read(map_file)
  if(trim(mask_file).eq."" .and. trim(imask_file).eq."" .and. trim(polmask_file).eq."") mask_file = "NONE"
  call sto%init(hot, peak_name, orient_name, nmaps = hgm%nmaps)
  sto%abs_threshold = hot .and. cold
  sto%angzero = .false.  
  sto%addpi = .true.
  sto%norm_power = 0.d0  !!these default settings can be changed when doing stacking
  select case(trim(coop_str_numUpperalpha(peak_name)))
     !!do nothing
  case("T", "I", "E", "B", "ZETA", "Z", "RANDOM", "SADDLE", "COL")
     if(hot)then
        sto%I_lower_nu = threshold
     else
        sto%I_upper_nu = -threshold
     endif
  case("P", "PT", "PZETA", "PZ")
     if(hot)then
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
  if(trim(mask_file).ne."NONE")then
     if(trim(mask_file).ne."")then
        call mask%read(mask_file, nmaps_wanted = 1)
     elseif( sto%mask_int .and. .not. sto%mask_pol)then
        call mask%read(imask_file, nmaps_wanted = 1)
     elseif( sto%mask_pol .and. .not. sto%mask_int)then
        call mask%read(polmask_file, nmaps_wanted = 1)
     else
        stop "For unknown class of peaks you need to specify the mask file explicitly"
     endif
     if(mask%nside .ne. hgm%nside) stop "mask and map must have the same nside"
     call hgm%apply_mask(mask = mask, remove_monopole = remove_mono)
     if(nside_scan .eq. 0)then
        call hgm%get_peaks(sto, mask = mask)
     else
        call hgm%get_peaks(sto, mask = mask, nside_scan = nside_scan)
     endif
  else
     if(nside_scan .eq. 0)then
        call hgm%get_peaks(sto)
     else
        call hgm%get_peaks(sto, nside_scan = nside_scan)
     endif
  endif
  print*, "find "//COOP_STR_OF(sto%peak_pix%n)//" spots"
  print*, "writing to "//trim(output)
  call sto%export(trim(output))
#else
  print*, "You need to install healpix"
#endif  
end program Exp_spots
