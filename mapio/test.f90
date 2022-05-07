program udg
  use coop_wrapper_utils
  use coop_healpix_mod
#ifdef HAS_HEALPIX  
  use udgrade_nr
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
#endif  
  implicit none
#include "constants.h"
#ifdef HAS_HEALPIX  
  character(LEN=1024)::f_map_in
  real*8 fac
  type(coop_healpix_maps)::map_in
  if(iargc().lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./Test -inp mapfile -fac 1.e6"
     stop
  endif
  call coop_get_command_line_argument(key = "inp", arg  = f_map_in)
  call coop_get_command_line_argument(key = "fac", arg = fac)
  call map_in%read(f_map_in)
  map_in%map = map_in%map * fac
  call map_in%write(f_map_in)
#else
  stop "you need to install healpix"
#endif  
end program udg
