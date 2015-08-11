program Test
  use coop_wrapper_typedef
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_ode)::ode
  call ode%init(n = 2)
  call ode%set_initial_conditions(xini = 0.d0, yini = (/ 1.d0, 0.d0 /))
  call ode%evolve(fcn, coop_pio2*3.d0)
  print*, ode%y
contains
  subroutine fcn(n, x, y, yp)
    COOP_INT::n
    COOP_REAL::x, y(n), yp(n)
    yp(1) =- (y(1)-cos(x))*100.d0
    yp(2) =- (y(2)+sin(x))*100.d0    
  end subroutine fcn
  
end program Test
