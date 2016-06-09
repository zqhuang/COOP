program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  COOP_INT,parameter::n = 256
  type(coop_function)::f
  COOP_REAL,dimension(n)::x, y
end program Test  
