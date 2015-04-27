program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_INT,parameter::n=10
  COOP_REAL::x(n) = (/ -1., 2., 3.2, 4.2, 5.1, 1.2, -2.3, -1.8, 4.44, 3.23 /)
  COOP_INT::minlocs(4), maxlocs(4)
  call coop_find_minlocs(x, minlocs)
  call coop_find_maxlocs(x, maxlocs)
  print*,minlocs
  print*, x(minlocs)
  print*, maxlocs
  print*, x(maxlocs)
  
end program Test
