program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_REAL::xmin, xmax, x
  x = 10.3d0
  xmin = -coop_pi
  xmax = coop_pi
  call coop_periodic_select_range(x, xmin, xmax)
  print*, x, x + coop_2pi*2


end program Test  
