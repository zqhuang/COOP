program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_REAL:: fnu
  type(coop_arguments)::args
  COOP_REAL,parameter::sigma0 = 1.d0, sigma1 = 1.d0, sigma2 = 2.d0
  COOP_INT,parameter::dim = 1
  call coop_gaussian_npeak_set_args(args, dim, sigma0, sigma1, sigma2)
  fnu = coop_integrate(coop_gaussian_npeak_differential, -6.d0, 6.d0, -6.d0, 0.d0, args, 1.d-5)
  print*, fnu, sigma2/sigma1/coop_2pi ! (sigma2/sigma1)**2/coop_pi/8.d0/coop_sqrt3
end program Test
