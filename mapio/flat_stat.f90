program test
  use coop_wrapper_utils
  use coop_fitswrap_mod
  use coop_healpix_mod
  use coop_sphere_mod
  implicit none
#include "constants.h"
  type(coop_fits_image_cea)::map
  COOP_STRING::fname
  call coop_MPI_Init()
  fname = trim(coop_InputArgs(1))
  if(trim(fname).eq."")then
     write(*,*) "Syntax: ./FStat FILENAME"
     stop
  endif
  call map%open(fname)
  call map%simple_stat()


end program test
