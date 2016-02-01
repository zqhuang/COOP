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
     write(*,*) "./FQU2EB QUMAP EBMAP"
     stop
  endif
  qumap = trim(coop_InputArgs(1))
  ebmap = trim(coop_InputArgs(2))
  call maps%read(qumap)
  if(.not. all(maps%spin .eq. 2)) stop "QU map must have spin 2"
  call maps%qu2eb()
  call maps%write(trim(ebmap)//"_EB.fsm")
  call maps%write(trim(ebmap)//"_E.fsm", write_image = .false., index_list = (/ 1 /) )
  call maps%write(trim(ebmap)//"_B.fsm", write_image = .false., index_list = (/ 2 /) )
end program test
