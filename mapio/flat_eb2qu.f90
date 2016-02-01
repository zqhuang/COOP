program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  COOP_REAL fwhm_arcmin
  COOP_STRING::qumap, ebmap
  type(coop_flatsky_maps)::maps
  if(iargc().lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./FQU2EB EBMAP QUMAP"
     stop
  endif
  ebmap = trim(coop_InputArgs(1))
  qumap = trim(coop_InputArgs(2))
  call maps%read(ebmap)
  call maps%eb2qu()
  call maps%write(qumap)
end program test
