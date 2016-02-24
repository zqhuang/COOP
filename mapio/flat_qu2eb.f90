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
  type(coop_fits_image_cea)::emap
  logical overwrite
  COOP_REAL::threshold, lambda
  if(iargc().lt.2)then
     write(*,*) "Syntax:"
     write(*,*) "./FQU2EB QUMAP EBMAP"
     stop
  endif
  qumap = trim(coop_InputArgs(1))
  ebmap = trim(coop_InputArgs(2))
  if(iargc() .gt. 2)then
     overwrite = trim( coop_inputargs(3)) == "T"
  else
     overwrite = .false.
  endif
  if(coop_file_exists(trim(ebmap)//"_EB.fsm") .and. .not. overwrite)then
     write(*,*) "the output file "//trim(ebmap)//"_EB.fsm already exists"
  else
     call maps%read(qumap)
     if(.not. all(maps%spin .eq. 2) .or. maps%nmaps .ne. 2) stop "QU map must have spin 2"
     call maps%qu2eb()
     !!try to do a global rotation do remove EB systematics
     !! E-> E+ lambda B , B -> B - lambda E
     !! so EB -> EB + lambda(B^2-E^2)
!!$     lambda = sum(maps%map(1)%image*maps%map(2)%image)/sum(maps%map(1)%image**2 - maps%map(2)%image**2)
!!$     write(*,*) "EB global rotation angle lambda = ", lambda
!!$     emap = maps%map(1)
!!$     maps%map(1)%image = emap%image + lambda*maps%map(2)%image
!!$     maps%map(2)%image = maps%map(2)%image - lambda * emap%image
     call maps%write(trim(ebmap)//"_EB.fsm")
     call maps%write(trim(ebmap)//"_E.fsm", write_image = .false., index_list = (/ 1 /) )
     call maps%write(trim(ebmap)//"_B.fsm", write_image = .false., index_list = (/ 2 /) )
  endif
end program test
