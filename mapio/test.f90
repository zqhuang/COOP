program shells
  use coop_fitswrap_mod
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT::i
  type(coop_fits_image_cea)::test
  call test%open("trial.fits")
end program shells
