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
  COOP_STRING::input, input2, output, genre, junk
  COOP_INT:: l, il, i
  type(coop_cls)::cls, cls2
  if(iargc() .lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./ClFormat -inp ... [-inp2 ...] -out ..."
     stop
  endif
  call coop_get_command_line_argument(key = "inp", arg = input)
  call coop_get_command_line_argument(key = "inp2", arg = input2, default = "")
  call coop_get_command_line_argument(key = "out", arg = output)
  call cls%load(input)
  if(trim(input2).ne."")then
     call cls2%load(input2)
     if(trim(cls%genre) .ne. trim(cls2%genre))then
        write(*,*) trim(cls%genre), " ", trim(cls2%genre)
        write(*,*) "inp and inp2 have different genre"
     endif
     call cls2%filter(lmin = cls%lmin, lmax = cls%lmax)
     cls%cls = cls%cls + cls2%cls
  endif
  call system('rm -f '//trim(output))
  call cls%dump(output)
end program txt2fits
