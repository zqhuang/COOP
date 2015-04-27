program test
  use coop_wrapper_firstorder
  use coop_forecast_mod
  implicit none
#include "constants.h"
  COOP_INT,parameter::n=399
  COOP_REAL,parameter :: alpha_min=0.8005d0, alpha_max = 1.1985d0
  type(coop_file)::fp
  COOP_REAL::alpha(n), prob(n)
  COOP_INT::i

end program test
