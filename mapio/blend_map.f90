program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::hp1, hp2
  COOP_INT::l_truncate
  COOP_STRING::map1, map2, output
  if(iargc() .lt. 4)then
     write(*,*) "Syntax:"
     write(*,*) "./Blend map1 map2 l_truncate output"
     write(*,*) "Example:"
     write(*,*) "./Blend planck_i.fits act_i.fits 500 planck_act_500_blend.fits"
     stop
  endif
  map1 = trim(coop_InputArgs(1))
  map2 = trim(coop_InputArgs(2))
  call coop_get_input(3, l_truncate)
  output = trim(coop_InputArgs(4))
  call hp1%read(map1)
  call hp2%read(map2)
  if(hp1%nside .ne. hp2%nside .or. hp1%nmaps .ne. hp2%nmaps) stop "cannot blend maps with different nside/nmaps"
  call hp1%map2alm()
  call hp2%map2alm()
  if(l_truncate .le. hp1%lmax)then
     hp1%alm(l_truncate:hp1%lmax,:,:) = hp2%alm(l_truncate:hp1%lmax,:,:)
     call hp1%alm2map()
  endif
  call hp1%write(trim(output))
end program test
