program test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"

  type(coop_function) cf
  type(coop_arguments) args
  args  = coop_arguments( r = (/ 1.d0/2.d0 /) )
  cf = coop_function(func, 0.1d0, 10.d0, xlog=.false., ylog=.false., args = args)
  print*, cf%derivative(2.d0)
contains
    
  function func(x, arg)
    COOP_REAL x, func
    type(coop_arguments) arg
    func = x**2
    func = func*arg%r(1)
  end function func
  
end program test
