program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  type(coop_ode)::ode
  call ode%init( n = 2, 
end program Test  
