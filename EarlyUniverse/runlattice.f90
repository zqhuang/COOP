program Test
  use coop_wrapper_utils
  use coop_lattice_mod
  implicit none
#include "constants.h"
#include "lattice.h"
  type(coop_lattice)::box
  call box%alloc(dim_f = 2, n = 512, nc = 128)
  box%pi_f(1, :, :, :) = 1.d0
  box%pi_f(2, :, :, :) = 2.d0
  print*, box%cs, box%nc
  call box%coarse_grain_pi()
  print*, box%pi_cf(1, 3, 2, 3), box%pi_cf(2, 3, 4,2)
end program Test
