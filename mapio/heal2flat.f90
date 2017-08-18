program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_STRING::input, output, sky
  COOP_REAL::ra_min, ra_max, dec_min, dec_max, pixsize
  type(coop_flatsky_maps)::mapout
  type(coop_healpix_maps)::mapin
  COOP_INT::i
  stop "not done yet"
  if(iargc().lt. 6)then
     write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++"
     write(*,*) "Convert healpix maps to RA-DEC flatsky maps"
     write(*,*) "Syntax:"
     write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++"     
     write(*,*) "./Heal2F -inp input.fits -out output.fsm -sky map1.fits"
     write(*,*) "convert the maps in file inputs.fits to a group of RA-DEC maps output.fsm (the postfix is for group of flat sky maps), the sky cut and map resolution are defined by map1.fits. "
     
     write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++"
     write(*,*) "-sky ... can be replaced with ranges of RA-DEC and pixsize, for example:"     
     write(*,*) "./Heal2F -inp input.fits -out output.fsm -RAmin 20 -RAmax 30 -DECmin -20 -DECmax 10 -pixsize 0.1666667"
     write(*,*) "convert the maps in file inputs.fits to a group of RA-DEC maps output.fsm; the sky cut is defined by 20<RA<30, -20<DEC<10; the resolution is 0.1666667 degree"
     write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++"
     stop
  endif
  call coop_get_command_line_argument(key = "inp", arg=input)
  call coop_get_command_line_argument(key = "out", arg=output)
  call coop_get_command_line_argument(key = "sky", arg=sky, default="")
  
  call mapin%read(input)
  if(trim(sky).eq."")then
     call coop_get_command_line_argument(key = "RAmin", arg=ra_min)
     call coop_get_command_line_argument(key = "RAmax", arg=ra_max)
     call coop_get_command_line_argument(key = "DECmax", arg=dec_max)
     call coop_get_command_line_argument(key = "DECmin", arg=dec_min)
     call coop_get_command_line_argument(key = "pixsize", arg=pixsize)
     call mapout%init(nmaps = mapin%nmaps, ra_min = ra_min, ra_max =ra_max, dec_min = dec_min, dec_max = dec_max, pixsize= pixsize, spin = mapin%spin)
  else
     call mapout%read_from_one(filename = sky, nmaps =mapin%nmaps )
  endif
  mapout%spin = mapin%spin
  mapout%fields = mapin%fields
  mapout%units = mapin%units
  do i=1, mapin%nmaps
     call mapout%map(i)%from_healpix(mapin, i)
  enddo
  call mapout%write(output, write_image = .true.)
end program test
