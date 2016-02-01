program Exp_spots
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_stacking_mod
  use coop_fitswrap_mod
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX
  logical::do_max 
  COOP_STRING::map_file, peak_name, orient_name, output
  type(coop_stacking_options)::sto
  type(coop_flatsky_maps)::fm
  COOP_REAL::threshold
  COOP_INT::i

  if(iargc().le.0)then
     write(*,"(A)") "----------------------------------------------------------"     
     write(*,"(A)") "Syntax:"
     write(*,"(A)") "./FGetPeaks -map MAP_FILE"
     write(*,"(A)") "----------------------------------------------------------"
     write(*,"(A)") "other options are:"     
     write(*,"(A)") "-out OUTPUT_FILE"
     write(*,"(A)") "-hot [T|F]"
     write(*,"(A)") "-peak [SADDLE|COL|RANDOM|T|E|B|P_T|P|\zeta|P_\zeta]"
     write(*,"(A)") "-orient [RANDOM|(Q_T, U_T)|(Q, U)|(Q_{\nabla^2T},U_{\nabla^2T})]"
     write(*,"(A)") "-nu THRESHOLD_VALUE"
     write(*,"(A)") "----------------------------------------------------------"     
     stop
  endif
  call coop_get_command_line_argument(key = 'map', arg = map_file)
  call coop_get_command_line_argument(key = 'hot', arg = do_max, default =.true.)

  call coop_get_command_line_argument(key = 'peak', arg = peak_name, default = 'T')
  call coop_get_command_line_argument(key = 'orient', arg = orient_name, default = 'RANDOM')
  call coop_get_command_line_argument(key = 'nu', arg = threshold, default = 0.d0)
  if((index(orient_name, coop_backslash).ne.0 .or. index(orient_name, "_") .ne. 0 .or. index(orient_name, "^").ne.0) .and. index(orient_name, "$").eq.0)then
     orient_name = "$"//trim(adjustl(orient_name))//"$"
  endif
  if(do_max)then
     call coop_get_command_line_argument(key = 'out', arg = output, default = "peaks/"//trim(coop_file_name_of(map_file, want_ext = .false.))//"_"//trim(coop_str_numalpha(peak_name))//"_"//trim(coop_str_numalpha(orient_name))//"_hot.dat")
  else
     call coop_get_command_line_argument(key = 'out', arg = output, default = "peaks/"//trim(coop_file_name_of(map_file, want_ext = .false.))//"_"//trim(coop_str_numalpha(peak_name))//"_"//trim(coop_str_numalpha(orient_name))//"_cold.dat")
  endif
  coop_healpix_mask_tol = 0.
  call fm%read(map_file)
  call sto%init(do_max, peak_name, orient_name, nmaps = fm%nmaps)
  sto%angzero = .false.  
  select case(trim(coop_str_numUpperalpha(peak_name)))
     !!do nothing
  case("T", "I", "E", "B", "ZETA", "Z", "RANDOM", "SADDLE", "COL")
     if(do_max)then
        sto%I_lower_nu = threshold
     else
        sto%I_upper_nu = -threshold
     endif
  case("P", "PT", "PZETA", "PZ")
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
  call fm%get_peaks(sto)
  print*, "find "//COOP_STR_OF(sto%peak_pix%n)//" spots"
  print*, "writing to "//trim(output)
  call sto%export(trim(output))
#else
  print*, "You need to install healpix"
#endif  
end program Exp_spots
