program txt2fits
  use coop_wrapper_firstorder
  use coop_zeta3d_mod
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  use coop_fitsio_mod
  use head_fits
  use fitstools
  use pix_tools
  use alm_tools
  implicit none
#include "constants.h"
  COOP_STRING::input, output, genre, junk
  COOP_INT:: l, il, i
  type(coop_cls)::cls
  if(iargc() .lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./ClFormat -inp ... -out ..."
     stop
  endif
  call coop_get_command_line_argument(key = "inp", arg = input)
  call coop_get_command_line_argument(key = "out", arg = output)
  call cls%load(input)
  call system('rm -f '//trim(output))
  call cls%dump(output)
end program txt2fits
