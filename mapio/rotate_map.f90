program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::hp
  COOP_STRING::fin, fout, coorin, coorout
  if(iargc().lt.4)then
     write(*,*) "Syntax:"
     write(*,*) "./Rotate input_file output_file  input_coor output_coor"
     write(*,*) "Example:"
     write(*,*) "./Rotate planck_i.fits planck_C_i.fits G C"
     stop
  endif
  fin = trim(coop_InputArgs(1))
  fout = trim(coop_InputArgs(2))
  coorin = trim(coop_InputArgs(3))
  coorout = trim(coop_InputArgs(4))
  call hp%read(fin)
  call hp%rotate_coor(trim(coorin), trim(coorout))
  call hp%write(fout)
end program test
