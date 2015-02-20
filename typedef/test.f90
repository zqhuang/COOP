program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  type(coop_function)::pol
  call pol%init_polynomial( (/ 0.5d0, 3.d0, 2.d0 /) )
  print*, pol%f2, pol%n
  print*, pol%eval(3.d0), pol%derivative(3.d0), pol%derivative2(3.d0)

end program Test
