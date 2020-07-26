program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_sphere_mod
  use coop_healpix_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::hp1, hp2, mask
  COOP_INT::ind
  COOP_REAL::summask
  COOP_STRING::map1, map2,fmask
  if(iargc() .lt. 3)then
     write(*,*) "Syntax:"
     write(*,*) "./Corr map1 map2 mask index_map"
     write(*,*) "Example:"
     write(*,*) "./Corr planck_i.fits act_i.fits mask.fits 1"
     stop
  endif
  map1 = trim(coop_InputArgs(1))
  map2 = trim(coop_InputArgs(2))
  fmask = trim(coop_InputArgs(3))
  if(iargc().gt.3)then
     call coop_get_input(4, ind)
  else
     ind = 1
  endif
  call hp1%read(map1)
  call hp2%read(map2)
  call mask%read(fmask, nmaps_wanted = 1)
  write(*,"(A4, F10.4)") "r = ", sum(dble(hp1%map(:,ind)*hp2%map(:,ind)*mask%map(:,1))) /sqrt( sum(dble(hp1%map(:,ind)**2 *mask%map(:,1)))  *  sum(dble(hp2%map(:,ind)**2*mask%map(:,1))) )
end program test
