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
  type(coop_file)::fp
  COOP_INT::binsize
  if(iargc() .lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./ClFormat -inp ... [-inp2 ... -binsize ... ] -out ..."
     write(*,*) "Example:"
     write(*,*) "./ClFormat -inp rawcls.fits -binsize 10 -out smoothedcls.dat (SMOOTH the cls in rawcls.fits with gaussian window (width~10) and write the smoothed cls to smoothedcls.dat)"
     stop
  endif
  call coop_get_command_line_argument(key = "inp", arg = input)
  call coop_get_command_line_argument(key = "inp2", arg = input2, default = "")
  call coop_get_command_line_argument(key = "out", arg = output)
  call coop_get_command_line_argument(key = "binsize", arg = binsize, default = 0)
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
  if(binsize .gt. 1)then
     call cls%smooth(binsize)
  endif
  call cls%dump(output)
end program txt2fits
