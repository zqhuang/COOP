program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  type(coop_ode)::ode
  call ode%init(n = 3)
  call ode%set_initial_conditions( xini = 0.d0, yini = (/ 0.d0, 1.d0, 0.d0 /) )
  call ode%evolve(get_yprime, coop_pi*4.d0)
  print*, ode%y
contains

  subroutine get_yprime(n, x, y, yp)
    COOP_INT n
    COOP_REAL x, y(0:n-1), yp(0:n-1)
    yp(0) = 0.d0
    yp(1) = y(2) + y(0)*100.d0
    yp(2) = - y(1) 
    if(y(0).ne.0.d0)then
       call coop_return_error("get_yprime", "nonzero y(0)", "stop")
    endif
  end subroutine get_yprime

end program Test
