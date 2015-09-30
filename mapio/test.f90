program test
  use coop_wrapper_utils
  use coop_healpix_mod
  use coop_fitswrap_mod
  use coop_stacking_mod
  implicit none
#include "constants.h"
  type(coop_healpix_maps)::map
  type(coop_healpix_disc) disc
  COOP_INT::pix
  COOP_REAL
  call coop_MPI_init()
  call map%read(
  call coop_MPI_Finalize()
end program test
