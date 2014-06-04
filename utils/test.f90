program plot
  use coop_wrapper_utils
  use coop_ode_mod
  implicit none
#include "constants.h"
  type(coop_arguments) args
  type(coop_ode) co
  integer,parameter::n = 2, nsteps = 100
  integer i
  args = coop_arguments( r = (/ coop_pi**2 /) )
  call co%init(n =n, method = COOP_ODE_GL6)
  call co%set_arguments(args)
  call co%set_initial_conditions(xini = 0.d0, yini = (/ 1.d0, 0.d0 /) )
  print*, co%y
  do i= 1, nsteps
     call co%evolve(fcn = get_yprime, xend = 9.5d0*i/nsteps)
  enddo
  print*, co%y
contains

  subroutine get_yprime(n, x, y, yp, args)
    COOP_INT n
    COOP_REAL x, y(n), yp(n)
    type(coop_arguments) args
    yp(1) = y(2)
    yp(2) = - args%r(1) * y(1)
  end subroutine get_yprime

end program plot
