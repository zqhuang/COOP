program Test
  use coop_wrapper_typedef
  implicit none
#include "constants.h"
  type(coop_function)::f
  COOP_REAL::x
  call f%init_rational(c_up = (/ 2.d0, 2.d0 /), alpha_up = (/ 2.d0, 0.d0 /), c_down = (/ 1.d0, 1.d0 /), alpha_down = (/ 1.d0, 0.d0 /))
  read(*,*) x
  print*, x,  f%eval(x), f%derivative(x), f%derivative2(x)
end program Test  
