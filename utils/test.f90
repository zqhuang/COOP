program Test
  use coop_wrapper_utils
  implicit none
#include "constants.h"
  COOP_INT, parameter::n = 200
  COOP_INT i
  COOP_REAL x(n), y(n), a, b
  type(coop_asy)::fig
  a = 0.d0
  b = 8.d0
  call coop_set_uniform(n, x, a, b)
  y = sin(x)
  call fig%open("test.txt")
  call fig%init()
  call fig%curve(x, y)
  call fig%add_legend("Good legend", "red", "solid", 5.)
  call fig%close()
end program Test
