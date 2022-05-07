program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_STRING::input, output, sky
  COOP_INT::nside
  type(coop_healpix_maps)::mapout, mask
  type(coop_flatsky_maps)::mapin
  COOP_INT::i
  stop "not done yet"
  if(iargc().lt. 8)then
     write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++"
     write(*,*) "Convert healpix maps to RA-DEC flatsky maps"
     write(*,*) "Syntax:"
     write(*,*) "./F2Heal -inp input.fsm -out output.fits -sky mask.fits -nside 512"
     write(*,*) "convert the maps in file inputs.fsm (fsm file points to a group of RA-DEC maps) to Healpix maps, with NSIDE=512. Save the healpix maps to output.fits; save the sky cut to mask.fits"
     write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++"          
     stop
  endif
  call coop_get_command_line_argument(key = "inp", arg=input)
  call coop_get_command_line_argument(key = "out", arg=output)
  call coop_get_command_line_argument(key = "sky", arg=sky)  
  call coop_get_command_line_argument(key = "nside", arg=nside)  
  call mapin%read(input)
  call mapin%To_Healpix(nside, mapout, mask)
  call mapout%write(output)
  call mask%write(sky)
end program test
