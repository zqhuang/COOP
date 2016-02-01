program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_STRING::fimap, fqmap, fumap, fmask, orient_name, peak_name, output
  type(coop_fits_image_cea)::imap, qmap, umap, mask
  logical::do_max, hasi, hasqu, hasmask
  COOP_INT::nmaps
  COOP_REAL::nu
  call coop_get_command_line_argument(key = 'imap', arg=fimap, default="")
  call coop_get_command_line_argument(key = 'qmap', arg=fqmap, default="")
  call coop_get_command_line_argument(key = 'umap', arg=fqmap, default="")
  hasi = trim(fimap).ne.""
  hasqu =trim(fqmap) .ne. "" .and trim(fump).ne.""
  nmaps = 0
  if(hasi) nmaps = nmaps + 1
  if(hasqu) nmaps = nmaps+2
  call coop_get_command_line_argument(key = 'mask', arg=fmask, default="")
  hasmask = trim(fmask).ne.""
  call coop_get_command_line_argument(key = 'nu', arg = nu, default = 0.d0)
  call coop_get_command_line_argument(key = 'orient', arg = orient_name, default="RANDOM")
  if((index(orient_name, coop_backslash).ne.0 .or. index(orient_name, "_") .ne. 0 .or. index(orient_name, "^").ne.0) .and. index(orient_name, "$").eq.0)then
     orient_name = "$"//trim(adjustl(orient_name))//"$"
  endif
  call coop_get_command_line_argument(key = 'peak', arg = peak_name, default="T")
  call coop_get_command_line_argument(key = 'hot', arg = do_max, default=.true.)
  call coop_get_command_line_argument(key = 'out', arg = output)
  call sto%init(do_max, peak_name, orient_name, nmaps = nmaps)
  sto%angzero = .false.
  select case(trim(coop_str_numUpperalpha(peak_name)))
     !!do nothing
  case("T", "I", "E", "B", "ZETA", "Z", "RANDOM", "SADDLE", "COL")
     if(do_max)then
        sto%I_lower_nu = nu
     else
        sto%I_upper_nu = -nu
     endif
  case("P", "PT", "PZETA", "PZ")
     if(do_max)then
        sto%P_lower_nu = max(nu, 0.d0)
     else
        if(nu .le. 0.d0) stop "For Pmin you must use positive nu"
        sto%P_upper_nu = nu
     endif
  case default
     if(abs(nu).lt. coop_stacking_max_nu)then
        write(*,*) "cannot apply threshold for unknown peak name: "//trim(peak_name)
        stop
     endif
  end select

  if(hasmask)then
     call mask%open(fmask)

end program test
