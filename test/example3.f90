program test
  use coop_wrapper
  implicit none
#include "constants.h"
  COOP_INT, parameter:: n = 128
  COOP_REAL:: x(n), y(n)
  type(coop_asy)::figure
  call coop_set_uniform(n, x, 0.d0, coop_2pi)
  y = sin(x)
  call figure%open("testfig.txt")
  call figure%init(xlabel = "$x$", ylabel = "$\sin x$")
  call figure%plot(x, y)
  call figure%close()
end program test
